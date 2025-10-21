# ============================================================
# MAPAS DE AUTOCORRELACIÓN ESPACIAL CON LÍMITES DEPARTAMENTALES
# Encuesta Nacional Agropecuaria 2024 - Perú
# ============================================================

# PASO 1: INSTALAR Y CARGAR LIBRERÍAS
# ============================================================

paquetes <- c("haven", "dplyr", "sf", "spdep", "ggplot2", "rnaturalearth", 
              "rnaturalearthdata", "viridis", "scales", "tidyr", "geodata")

paquetes_faltantes <- paquetes[!(paquetes %in% installed.packages()[,"Package"])]
if(length(paquetes_faltantes) > 0) {
  install.packages(paquetes_faltantes)
}

library(haven)
library(dplyr)
library(sf)
library(spdep)
library(ggplot2)
library(viridis)
library(scales)

cat("✓ Librerías cargadas\n\n")

# ============================================================
# PASO 2: DESCARGAR MAPA DE PERÚ (AUTOMÁTICO)
# ============================================================

cat("Descargando mapa de Perú...\n")

# Opción 1: Usar rnaturalearth (más fácil pero menos detallado)
library(rnaturalearth)
library(rnaturalearthdata)

# Descargar países de Sudamérica
sudamerica <- ne_countries(continent = "south america", returnclass = "sf", scale = "medium")

# Filtrar solo Perú
peru_pais <- sudamerica %>% filter(name == "Peru")

# Descargar divisiones administrativas de Perú (nivel 1 = departamentos)
# Esto puede tardar un poco la primera vez
library(geodata)

# Descargar shapefile de departamentos del Perú
peru_dept_gadm <- gadm(country = "PER", level = 1, path = tempdir())
peru_dept <- st_as_sf(peru_dept_gadm)

# Limpiar nombres para hacer match con tus datos
peru_dept <- peru_dept %>%
  mutate(
    NOMBREDD_CLEAN = toupper(NAME_1),
    NOMBREDD_CLEAN = case_when(
      NOMBREDD_CLEAN == "LIMA METROPOLITANA" ~ "LIMA",
      NOMBREDD_CLEAN == "CALLAO" ~ "LIMA",
      TRUE ~ NOMBREDD_CLEAN
    )
  )

cat("✓ Mapa de Perú descargado\n")
cat("  Departamentos encontrados:", nrow(peru_dept), "\n\n")

# ============================================================
# PASO 3: CARGAR TUS DATOS
# ============================================================

# AJUSTA ESTAS RUTAS
ruta_cultivos <- "C:/Users/ROSMERY/Downloads/2024/2024 - DESCOMPRIMIDO/973-Modulo1893/01_CAP100A_04.sav"
ruta_caratula <- "C:/Users/ROSMERY/Downloads/2024/2024 - DESCOMPRIMIDO/973-Modulo1893/CARATULA.sav"

cat("Cargando datos...\n")
cultivos <- read_sav(ruta_cultivos)
caratula <- read_sav(ruta_caratula)

# ============================================================
# PASO 4: PROCESAR DATOS (IGUAL QUE ANTES)
# ============================================================

superficie_dept <- cultivos %>%
  filter(!is.na(P117_SUP_ha), P117_SUP_ha > 0, !is.na(FACTOR_PRODUCTOR)) %>%
  group_by(NOMBREDD) %>%
  summarise(
    superficie_total_ha = sum(P117_SUP_ha * FACTOR_PRODUCTOR, na.rm = TRUE),
    n_productores_expan = sum(FACTOR_PRODUCTOR, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    sup_promedio = superficie_total_ha / n_productores_expan,
    log_superficie = log(superficie_total_ha + 1),
    NOMBREDD_CLEAN = toupper(trimws(NOMBREDD))
  )

coords_dept <- caratula %>%
  filter(!is.na(LATITUD), !is.na(LONGITUD)) %>%
  group_by(NOMBREDD) %>%
  summarise(
    lon = mean(LONGITUD, na.rm = TRUE),
    lat = mean(LATITUD, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(NOMBREDD_CLEAN = toupper(trimws(NOMBREDD)))

datos_completos <- coords_dept %>%
  inner_join(superficie_dept, by = "NOMBREDD_CLEAN") %>%
  filter(!is.na(lon), !is.na(lat), !is.na(superficie_total_ha))

# Crear objeto espacial
datos_sf <- st_as_sf(datos_completos, 
                     coords = c("lon", "lat"), 
                     crs = 4326,
                     remove = FALSE)

# ============================================================
# PASO 5: CALCULAR AUTOCORRELACIÓN ESPACIAL
# ============================================================

cat("Calculando autocorrelación espacial...\n")

# Matriz de vecindad
coords_matrix <- st_coordinates(datos_sf)
k_vecinos <- 4
knn_nb <- knn2nb(knearneigh(coords_matrix, k = k_vecinos))
knn_weights <- nb2listw(knn_nb, style = "W", zero.policy = TRUE)

# Moran's I Global
moran_test <- moran.test(datos_sf$superficie_total_ha, 
                         knn_weights, 
                         zero.policy = TRUE)

# LISA
lisa_result <- localmoran(datos_sf$superficie_total_ha, 
                          knn_weights, 
                          zero.policy = TRUE)

datos_sf$lisa_I <- lisa_result[, 1]
datos_sf$lisa_pvalue <- lisa_result[, 5]
datos_sf$z_superficie <- scale(datos_sf$superficie_total_ha)[,1]
lag_superficie <- lag.listw(knn_weights, datos_sf$superficie_total_ha)
datos_sf$z_lag <- scale(lag_superficie)[,1]

datos_sf$lisa_cluster <- case_when(
  datos_sf$lisa_pvalue >= 0.05 ~ "No significativo",
  datos_sf$z_superficie > 0 & datos_sf$z_lag > 0 ~ "High-High (HH)",
  datos_sf$z_superficie < 0 & datos_sf$z_lag < 0 ~ "Low-Low (LL)",
  datos_sf$z_superficie > 0 & datos_sf$z_lag < 0 ~ "High-Low (HL)",
  datos_sf$z_superficie < 0 & datos_sf$z_lag > 0 ~ "Low-High (LH)",
  TRUE ~ "No significativo"
)

# Gi*
gi_star <- localG(datos_sf$superficie_total_ha, 
                  knn_weights, 
                  zero.policy = TRUE)

datos_sf$gi_star <- as.numeric(gi_star)
datos_sf$hotspot <- case_when(
  gi_star > 2.58 ~ "Hotspot muy significativo (p<0.01)",
  gi_star > 1.96 ~ "Hotspot significativo (p<0.05)",
  gi_star < -2.58 ~ "Coldspot muy significativo (p<0.01)",
  gi_star < -1.96 ~ "Coldspot significativo (p<0.05)",
  TRUE ~ "No significativo"
)

cat("✓ Autocorrelación calculada\n\n")

# ============================================================
# PASO 6: UNIR DATOS CON MAPA DE DEPARTAMENTOS
# ============================================================

cat("Uniendo datos con mapa...\n")

# Unir shapefile con tus resultados
peru_con_datos <- peru_dept %>%
  left_join(
    datos_sf %>% st_drop_geometry(),
    by = "NOMBREDD_CLEAN"
  )

# Ver qué departamentos no hicieron match
sin_datos <- peru_con_datos %>% 
  filter(is.na(superficie_total_ha)) %>% 
  pull(NOMBREDD_CLEAN)

if(length(sin_datos) > 0) {
  cat("⚠️ Departamentos sin datos:\n")
  print(sin_datos)
  cat("\n")
}

cat("✓ Datos unidos con mapa\n\n")

# ============================================================
# PASO 7: CREAR MAPAS PROFESIONALES
# ============================================================

cat("Generando mapas...\n\n")

# Tema común para todos los mapas
tema_mapa <- theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(size = 9, hjust = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid = element_line(color = "gray90"),
    axis.text = element_text(size = 9)
  )

# ============================================================
# MAPA 1: SUPERFICIE AGRÍCOLA (COROPLÉTICO)
# ============================================================

cat("📍 Mapa 1: Superficie agrícola por departamento...\n")

mapa1 <- ggplot(peru_con_datos) +
  # Fondo del país
  geom_sf(aes(fill = superficie_total_ha / 1000), color = "white", size = 0.3) +
  # Escala de colores
  scale_fill_viridis_c(
    name = "Superficie\n(miles de ha)",
    option = "inferno",
    na.value = "gray90",
    labels = comma,
    trans = "log10"
  ) +
  # Etiquetas de departamentos (opcional)
  geom_sf_text(
    data = peru_con_datos %>% filter(!is.na(superficie_total_ha)),
    aes(label = NAME_1),
    size = 2.5,
    color = "white",
    fontface = "bold",
    check_overlap = TRUE
  ) +
  # Títulos
  labs(
    title = "Superficie Agrícola por Departamento",
    subtitle = "Encuesta Nacional Agropecuaria 2024",
    caption = "Fuente: ENA 2024 - INEI"
  ) +
  tema_mapa

print(mapa1)
ggsave("MAPA_01_superficie_agricola.png", mapa1, width = 10, height = 12, dpi = 300, bg = "white")
cat("✓ Guardado: MAPA_01_superficie_agricola.png\n\n")

# ============================================================
# MAPA 2: CLÚSTERES LISA
# ============================================================

cat("📍 Mapa 2: Clústeres LISA...\n")

mapa2 <- ggplot(peru_con_datos) +
  # Fondo del país
  geom_sf(aes(fill = lisa_cluster), color = "white", size = 0.3) +
  # Colores de clústeres
  scale_fill_manual(
    name = "Clúster LISA",
    values = c(
      "High-High (HH)" = "#d7191c",
      "Low-Low (LL)" = "#2c7bb6",
      "High-Low (HL)" = "#fdae61",
      "Low-High (LH)" = "#abd9e9",
      "No significativo" = "#f0f0f0"
    ),
    na.value = "gray90",
    breaks = c("High-High (HH)", "Low-Low (LL)", "High-Low (HL)", 
               "Low-High (LH)", "No significativo")
  ) +
  # Etiquetas solo para departamentos significativos
  geom_sf_text(
    data = peru_con_datos %>% filter(lisa_cluster != "No significativo", !is.na(lisa_cluster)),
    aes(label = NAME_1),
    size = 3,
    color = "black",
    fontface = "bold"
  ) +
  # Títulos
  labs(
    title = "Clústeres LISA - Análisis Local de Autocorrelación Espacial",
    subtitle = paste0("Superficie Agrícola | Moran's I = ", round(moran_test$estimate[1], 3), 
                      " (p = ", round(moran_test$p.value, 4), ")"),
    caption = "HH: Alta rodeada de alta | LL: Baja rodeada de baja\nHL/LH: Outliers espaciales"
  ) +
  tema_mapa

print(mapa2)
ggsave("MAPA_02_clusters_LISA.png", mapa2, width = 10, height = 12, dpi = 300, bg = "white")
cat("✓ Guardado: MAPA_02_clusters_LISA.png\n\n")

# ============================================================
# MAPA 3: HOTSPOTS Y COLDSPOTS (Gi*)
# ============================================================

cat("📍 Mapa 3: Hotspots y Coldspots...\n")

mapa3 <- ggplot(peru_con_datos) +
  # Fondo del país
  geom_sf(aes(fill = hotspot), color = "white", size = 0.3) +
  # Colores de hotspots
  scale_fill_manual(
    name = "Clasificación Gi*",
    values = c(
      "Hotspot muy significativo (p<0.01)" = "#b10026",
      "Hotspot significativo (p<0.05)" = "#fc4e2a",
      "Coldspot muy significativo (p<0.01)" = "#08519c",
      "Coldspot significativo (p<0.05)" = "#6baed6",
      "No significativo" = "#f0f0f0"
    ),
    na.value = "gray90"
  ) +
  # Etiquetas para hotspots/coldspots
  geom_sf_text(
    data = peru_con_datos %>% filter(hotspot != "No significativo", !is.na(hotspot)),
    aes(label = NAME_1),
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  # Títulos
  labs(
    title = "Hotspots y Coldspots - Getis-Ord Gi*",
    subtitle = "Concentraciones de Superficie Agrícola - ENA 2024",
    caption = "Hotspot: Concentración alta | Coldspot: Concentración baja"
  ) +
  tema_mapa

print(mapa3)
ggsave("MAPA_03_hotspots_gi.png", mapa3, width = 10, height = 12, dpi = 300, bg = "white")
cat("✓ Guardado: MAPA_03_hotspots_gi.png\n\n")

# ============================================================
# MAPA 4: VALORES GI* CONTINUOS
# ============================================================

cat("📍 Mapa 4: Valores Gi* continuos...\n")

mapa4 <- ggplot(peru_con_datos) +
  # Fondo del país
  geom_sf(aes(fill = gi_star), color = "white", size = 0.3) +
  # Escala divergente
  scale_fill_gradient2(
    name = "Estadístico\nGi*",
    low = "#2c7bb6",
    mid = "#ffffbf",
    high = "#d7191c",
    midpoint = 0,
    na.value = "gray90",
    limits = c(-3, 3),
    breaks = c(-2.58, -1.96, 0, 1.96, 2.58),
    labels = c("-2.58\n(p<0.01)", "-1.96\n(p<0.05)", "0", 
               "+1.96\n(p<0.05)", "+2.58\n(p<0.01)")
  ) +
  # Etiquetas
  geom_sf_text(
    data = peru_con_datos %>% filter(abs(gi_star) > 1.96, !is.na(gi_star)),
    aes(label = NAME_1),
    size = 2.5,
    color = "black",
    fontface = "bold"
  ) +
  # Títulos
  labs(
    title = "Valores del Estadístico Getis-Ord Gi*",
    subtitle = "Escala continua de concentración espacial",
    caption = "Valores > 1.96 o < -1.96 son estadísticamente significativos (p < 0.05)"
  ) +
  tema_mapa

print(mapa4)
ggsave("MAPA_04_gi_continuo.png", mapa4, width = 10, height = 12, dpi = 300, bg = "white")
cat("✓ Guardado: MAPA_04_gi_continuo.png\n\n")

# ============================================================
# MAPA 5: COMBINADO - SUPERFICIE + CLÚSTERES
# ============================================================

cat("📍 Mapa 5: Mapa combinado...\n")

mapa5 <- ggplot() +
  # Base: superficie agrícola
  geom_sf(
    data = peru_con_datos,
    aes(fill = superficie_total_ha / 1000),
    color = "white",
    size = 0.3
  ) +
  scale_fill_viridis_c(
    name = "Superficie\n(miles de ha)",
    option = "viridis",
    na.value = "gray90",
    trans = "log10",
    alpha = 0.7
  ) +
  # Overlay: Bordes de clústeres significativos
  geom_sf(
    data = peru_con_datos %>% filter(lisa_cluster == "High-High (HH)"),
    fill = NA,
    color = "#d7191c",
    size = 1.5,
    linetype = "solid"
  ) +
  geom_sf(
    data = peru_con_datos %>% filter(lisa_cluster == "Low-High (LH)"),
    fill = NA,
    color = "#abd9e9",
    size = 1.5,
    linetype = "dashed"
  ) +
  # Etiquetas
  geom_sf_text(
    data = peru_con_datos %>% filter(lisa_cluster %in% c("High-High (HH)", "Low-High (LH)")),
    aes(label = NAME_1),
    size = 3,
    color = "white",
    fontface = "bold"
  ) +
  # Títulos
  labs(
    title = "Superficie Agrícola y Clústeres Espaciales",
    subtitle = "Línea roja sólida: Clúster High-High | Línea celeste punteada: Low-High",
    caption = "Fuente: ENA 2024 - INEI"
  ) +
  tema_mapa

print(mapa5)
ggsave("MAPA_05_combinado.png", mapa5, width = 10, height = 12, dpi = 300, bg = "white")
cat("✓ Guardado: MAPA_05_combinado.png\n\n")

# ============================================================
# MAPA 6: PANEL DE 4 MAPAS
# ============================================================

cat("📍 Mapa 6: Panel comparativo...\n")

library(gridExtra)
library(grid)

# Crear versiones simplificadas para el panel
p1_panel <- ggplot(peru_con_datos) +
  geom_sf(aes(fill = superficie_total_ha / 1000), color = "white", size = 0.2) +
  scale_fill_viridis_c(name = "Sup. (miles ha)", option = "inferno", trans = "log10", na.value = "gray90") +
  labs(title = "A) Superficie Agrícola") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.text = element_blank(),
    panel.grid = element_blank()
  )

p2_panel <- ggplot(peru_con_datos) +
  geom_sf(aes(fill = lisa_cluster), color = "white", size = 0.2) +
  scale_fill_manual(
    name = "Clúster",
    values = c("High-High (HH)" = "#d7191c", "Low-High (LH)" = "#abd9e9", "No significativo" = "#f0f0f0"),
    na.value = "gray90"
  ) +
  labs(title = "B) Clústeres LISA") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.text = element_blank(),
    panel.grid = element_blank()
  )

p3_panel <- ggplot(peru_con_datos) +
  geom_sf(aes(fill = hotspot), color = "white", size = 0.2) +
  scale_fill_manual(
    name = "Gi*",
    values = c(
      "Hotspot muy significativo (p<0.01)" = "#b10026",
      "Hotspot significativo (p<0.05)" = "#fc4e2a",
      "No significativo" = "#f0f0f0"
    ),
    na.value = "gray90"
  ) +
  labs(title = "C) Hotspots Gi*") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.text = element_blank(),
    panel.grid = element_blank()
  )

p4_panel <- ggplot(peru_con_datos) +
  geom_sf(aes(fill = gi_star), color = "white", size = 0.2) +
  scale_fill_gradient2(
    name = "Gi* valor",
    low = "#2c7bb6", mid = "#ffffbf", high = "#d7191c",
    midpoint = 0, na.value = "gray90", limits = c(-3, 3)
  ) +
  labs(title = "D) Gi* Continuo") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    legend.key.height = unit(0.3, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.text = element_blank(),
    panel.grid = element_blank()
  )

# Combinar en panel
panel <- grid.arrange(
  p1_panel, p2_panel, p3_panel, p4_panel,
  ncol = 2,
  top = textGrob(
    "Análisis de Autocorrelación Espacial - Superficie Agrícola del Perú (ENA 2024)",
    gp = gpar(fontface = "bold", fontsize = 14)
  ),
  bottom = textGrob(
    "Fuente: Encuesta Nacional Agropecuaria 2024 - INEI",
    gp = gpar(fontsize = 9),
    hjust = 1,
    x = 0.95
  )
)

ggsave("MAPA_06_panel_comparativo.png", panel, width = 14, height = 16, dpi = 300, bg = "white")
cat("✓ Guardado: MAPA_06_panel_comparativo.png\n\n")

# ============================================================
# EXPORTAR SHAPEFILE CON RESULTADOS
# ============================================================

cat("Exportando shapefile con resultados...\n")

st_write(peru_con_datos, "PERU_autocorrelacion_espacial.shp", delete_dsn = TRUE)
cat("✓ Shapefile guardado: PERU_autocorrelacion_espacial.shp\n")
cat("  (Puedes abrirlo en QGIS, ArcGIS, etc.)\n\n")
#FIN
