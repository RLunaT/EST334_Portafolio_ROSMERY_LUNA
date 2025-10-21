# ============================================================================
# ANÁLISIS ESPACIAL DE PRODUCCIÓN AGRÍCOLA E IRRIGACIÓN EN PERÚ
# Preparación de Datos y Control de Calidad
# ============================================================================

# 1. CARGAR LIBRERÍAS NECESARIAS
# ============================================================================
library(haven)        # Ya cargado - para leer archivos SPSS
library(tidyverse)    # Manipulación de datos
library(sf)           # Datos espaciales (Simple Features)
library(sp)           # Clases espaciales
library(gstat)        # Geoestadística
library(spdep)        # Autocorrelación espacial
library(viridis)      # Paletas de colores
library(tmap)         # Mapas temáticos

library(haven)
CARATULA <- read_sav("C:/Users/ROSMERY/Desktop/PREGRADO/UNA-PUNO/FINESI/MIS SEMESTRES ACÁDEMICOS/10cmo/ESTADÍSTICA ESPACIAL/2024 - DESCOMPRIMIDO/973-Modulo1893/CARATULA.sav")

X03_CAP200AB <- read_sav("C:/Users/ROSMERY/Desktop/PREGRADO/UNA-PUNO/FINESI/MIS SEMESTRES ACÁDEMICOS/10cmo/ESTADÍSTICA ESPACIAL/2024 - DESCOMPRIMIDO/973-Modulo1895/03_CAP200AB.sav")


# 2. CONTROL DE CALIDAD Y FILTRADO DE DATOS
# ============================================================================

# 2.1 Filtrar observaciones completas y con coordenadas válidas
cat("Dimensiones originales CARATULA:", dim(CARATULA), "\n")
cat("Dimensiones originales CAP200AB:", dim(X03_CAP200AB), "\n")

# Filtrar CARATULA: solo entrevistas completas y coordenadas válidas
CARATULA_limpia <- CARATULA %>%
  filter(RESFIN == 1,  # 1 = Completa
         !is.na(LATITUD),
         !is.na(LONGITUD),
         LATITUD != 0,
         LONGITUD != 0) %>%
  # Crear identificador único de segmento-productor
  mutate(ID_SEGM_PROD = paste(NSEGM, ID_PROD, sep = "_"))

# Filtrar CAP200AB: solo entrevistas completas
CAP200AB_limpia <- X03_CAP200AB %>%
  filter(RESFIN == 1,  # 1 = Completa
         !is.na(P217_SUP_ha),    # Superficie cosechada válida
         !is.na(P219_EQUIV_KG),  # Producción válida
         P217_SUP_ha > 0,        # Superficie mayor a cero
         P219_EQUIV_KG > 0) %>%  # Producción mayor a cero
  mutate(ID_SEGM_PROD = paste(NSEGM, ID_PROD, sep = "_"))

cat("\nDimensiones después del filtrado:\n")
cat("CARATULA_limpia:", dim(CARATULA_limpia), "\n")
cat("CAP200AB_limpia:", dim(CAP200AB_limpia), "\n")

# 2.2 Unir ambas bases de datos
# Primero eliminamos TODAS las columnas que ya existen en CAP200AB_limpia
datos_completos <- CAP200AB_limpia %>%
  select(-c(REGION, ESTRATO, FACTOR_PRODUCTOR, NOMBREDD, NOMBREPV, NOMBREDI)) %>%  # Eliminar duplicadas
  left_join(CARATULA_limpia %>% 
              select(ID_SEGM_PROD, LATITUD, LONGITUD, FACTOR_PRODUCTOR,
                     REGION, ESTRATO, NOMBREDD, NOMBREPV, NOMBREDI),
            by = "ID_SEGM_PROD") %>%
  filter(!is.na(LATITUD), !is.na(LONGITUD))

cat("\nDatos unidos:", nrow(datos_completos), "observaciones\n")

# 3. AGREGAR DATOS A NIVEL DE SEGMENTO
# ============================================================================

# 3.1 Primero obtener las coordenadas y variables geográficas únicas por segmento
segmento_geo <- datos_completos %>%
  group_by(NSEGM) %>%
  summarise(
    LATITUD = first(LATITUD),
    LONGITUD = first(LONGITUD),
    REGION = first(REGION),
    ESTRATO = first(ESTRATO),
    NOMBREDD = first(NOMBREDD),
    NOMBREPV = first(NOMBREPV),
    NOMBREDI = first(NOMBREDI),
    .groups = "drop"
  )

# 3.2 Estadísticas por segmento
datos_segmento <- datos_completos %>%
  group_by(NSEGM) %>%
  summarise(
    # Producción total ponderada por factor de expansión
    prod_total_kg = sum(P219_EQUIV_KG * FACTOR_PRODUCTOR, na.rm = TRUE),
    
    # Superficie total cosechada ponderada
    sup_total_ha = sum(P217_SUP_ha * FACTOR_PRODUCTOR, na.rm = TRUE),
    
    # Intensidad de producción (kg/ha)
    intensidad_prod = ifelse(sup_total_ha > 0, 
                             prod_total_kg / sup_total_ha, 
                             NA),
    
    # Número de cultivos diferentes
    n_cultivos = n_distinct(P204_NOM),
    
    # Índice de diversidad de Shannon
    shannon_index = -sum((P217_SUP_ha / sum(P217_SUP_ha, na.rm = TRUE)) * 
                           log(P217_SUP_ha / sum(P217_SUP_ha, na.rm = TRUE) + 1e-10),
                         na.rm = TRUE),
    
    # Proporción con riego
    prop_riego = mean(P212 != 1, na.rm = TRUE),  # P212 != 1 significa que NO es lluvia
    
    # Fuente de agua predominante
    fuente_agua_pred = names(which.max(table(P212))),
    
    # Sistema de riego predominante (entre los que tienen riego)
    sistema_riego_pred = ifelse(any(!is.na(P213)), 
                                names(which.max(table(P213[!is.na(P213)]))), 
                                NA),
    
    # Factor de expansión promedio
    factor_exp = mean(FACTOR_PRODUCTOR, na.rm = TRUE),
    
    # Número de observaciones por segmento
    n_obs = n(),
    
    .groups = "drop"
  ) %>%
  # Unir con la información geográfica
  left_join(segmento_geo, by = "NSEGM")

cat("\nSegmentos únicos:", nrow(datos_segmento), "\n")

# 4. CREAR OBJETO ESPACIAL
# ============================================================================

# Convertir a objeto sf (Simple Features)
datos_sf <- st_as_sf(datos_segmento, 
                     coords = c("LONGITUD", "LATITUD"),
                     crs = 4326)  # WGS84

# Resumen de las variables
cat("\n=== RESUMEN DE VARIABLES CLAVE ===\n")
summary(datos_segmento[, c("prod_total_kg", "sup_total_ha", 
                           "intensidad_prod", "n_cultivos", 
                           "shannon_index", "prop_riego")])

# 5. ANÁLISIS EXPLORATORIO ESPACIAL
# ============================================================================

# 5.1 Distribución por región
cat("\n=== DISTRIBUCIÓN POR REGIÓN ===\n")
table(datos_segmento$REGION)

# 5.2 Estadísticas por región
stats_region <- datos_segmento %>%
  group_by(REGION) %>%
  summarise(
    n_segmentos = n(),
    prod_media_kg = mean(prod_total_kg, na.rm = TRUE),
    sup_media_ha = mean(sup_total_ha, na.rm = TRUE),
    intensidad_media = mean(intensidad_prod, na.rm = TRUE),
    prop_con_riego = mean(prop_riego > 0, na.rm = TRUE)
  )

print(stats_region)

# 5.3 Identificar valores atípicos
# Usar método de Tukey (Q1 - 1.5*IQR, Q3 + 1.5*IQR)
detect_outliers <- function(x) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  
  outliers <- x < (q1 - 1.5 * iqr) | x > (q3 + 1.5 * iqr)
  return(outliers)
}

datos_segmento <- datos_segmento %>%
  mutate(
    outlier_prod = detect_outliers(prod_total_kg),
    outlier_intensidad = detect_outliers(intensidad_prod)
  )

cat("\nValores atípicos detectados:\n")
cat("Producción:", sum(datos_segmento$outlier_prod, na.rm = TRUE), "\n")
cat("Intensidad:", sum(datos_segmento$outlier_intensidad, na.rm = TRUE), "\n")

# 6. ANÁLISIS DE IRRIGACIÓN
# ============================================================================

# 6.1 Categorizar fuente de agua
datos_irrigacion <- datos_completos %>%
  mutate(
    fuente_agua_cat = case_when(
      P212 == 1 ~ "Sin riego (Lluvia)",
      P212 == 2 ~ "Río",
      P212 == 3 ~ "Manantial/Puquio",
      P212 == 4 ~ "Pozo/Subterránea",
      P212 == 5 ~ "Reservorio/Laguna",
      P212 == 6 ~ "Represa",
      P212 == 7 ~ "Canal de riego",
      TRUE ~ "Otro"
    ),
    
    sistema_riego_cat = case_when(
      is.na(P213) ~ "Sin tecnología",
      P213 == 1 ~ "Exudación",
      P213 == 2 ~ "Goteo",
      P213 == 3 ~ "Microaspersión",
      P213 == 4 ~ "Aspersión",
      P213 == 5 ~ "Gravedad",
      P213 == 6 ~ "Inundación",
      P213 == 7 ~ "Otro",
      P213 == 8 ~ "Ninguno",
      TRUE ~ "Desconocido"
    )
  )

# Tabla de frecuencias
cat("\n=== FUENTE DE AGUA ===\n")
table(datos_irrigacion$fuente_agua_cat)

cat("\n=== SISTEMA DE RIEGO ===\n")
table(datos_irrigacion$sistema_riego_cat)

# 6.2 Agregar irrigación por segmento
irrigacion_segmento <- datos_irrigacion %>%
  group_by(NSEGM) %>%
  summarise(
    fuente_predominante = names(which.max(table(fuente_agua_cat))),
    sistema_predominante = names(which.max(table(sistema_riego_cat))),
    .groups = "drop"
  )

# Unir con datos espaciales
datos_segmento <- datos_segmento %>%
  left_join(irrigacion_segmento, by = "NSEGM")

# 7. GUARDAR DATOS PROCESADOS
# ============================================================================

# Guardar como archivo RDS (formato R nativo)
saveRDS(datos_segmento, "datos_segmento_procesados.rds")
saveRDS(datos_sf, "datos_espaciales_sf.rds")
saveRDS(datos_completos, "datos_completos_cultivos.rds")

# Guardar como CSV
write_csv(datos_segmento, "datos_segmento_procesados.csv")

# ============================================================================
# VISUALIZACIONES ESENCIALES - ANÁLISIS ESPACIAL PERÚ
# ============================================================================

# LIBRERÍAS
library(tidyverse)
library(sf)
library(ggspatial)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

# ============================================================================
# VISUALIZACIONES ESENCIALES - ANÁLISIS ESPACIAL PERÚ
# ============================================================================

# LIBRERÍAS
library(tidyverse)
library(sf)
library(ggspatial)
library(viridis)

# Si no las tienes:
# install.packages(c("sf", "ggspatial", "viridis"))

# ============================================================================
# PREPARAR SHAPEFILES DE PERÚ (GADM)
# ============================================================================

# Crear carpeta si no existe
if(!dir.exists("shapefiles")) dir.create("shapefiles")

# Descargar si no existe
if(!file.exists("shapefiles/gadm41_PER_shp.zip")) {
  cat("  Descargando shapefiles de Perú desde GADM...\n")
  download.file("https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_PER_shp.zip",
                "shapefiles/gadm41_PER_shp.zip", mode = "wb", quiet = TRUE)
  unzip("shapefiles/gadm41_PER_shp.zip", exdir = "shapefiles/")
  cat("  ✓ Descarga completada\n")
}

# Cargar shapefiles
peru_pais <- st_read("shapefiles/gadm41_PER_0.shp", quiet = TRUE)  # Límite nacional
peru_dep <- st_read("shapefiles/gadm41_PER_1.shp", quiet = TRUE)   # Departamentos

# Simplificar geometrías para mejor rendimiento
peru_pais_simple <- st_simplify(peru_pais, dTolerance = 0.01)
peru_dep_simple <- st_simplify(peru_dep, dTolerance = 0.01)

cat("  ✓ Shapefiles cargados:", nrow(peru_dep), "departamentos\n\n")

# ASEGURAR QUE DATOS_SF EXISTE
if(!exists("datos_sf")) {
  datos_sf <- st_as_sf(datos_segmento, 
                       coords = c("LONGITUD", "LATITUD"),
                       crs = 4326)
}

# ============================================================================
# MAPA 1: DISTRIBUCIÓN ESPACIAL DE SEGMENTOS
# ============================================================================

mapa_segmentos <- ggplot() +
  geom_sf(data = peru_dep_simple, fill = "white", color = "gray30", size = 0.5) +
  geom_sf(data = datos_sf, aes(color = as.factor(REGION)), size = 0.8, alpha = 0.6) +
  scale_color_manual(
    values = c("1" = "#E74C3C", "2" = "#3498DB", "3" = "#27AE60"),
    labels = c("1" = "Costa", "2" = "Sierra", "3" = "Selva"),
    name = "Región Natural"
  ) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    title = "Distribución Espacial de Segmentos Muestrales",
    subtitle = "Censo Nacional Agropecuario 2024 - Perú",
    caption = "Fuente: INEI"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    legend.position = "right"
  )

print(mapa_segmentos)


# ============================================================================
# MAPA 2: INTENSIDAD DE PRODUCCIÓN
# ============================================================================

datos_sf_intensidad <- datos_sf %>%
  filter(!is.na(intensidad_prod), 
         intensidad_prod > 0,
         intensidad_prod < quantile(intensidad_prod, 0.99, na.rm = TRUE))

mapa_intensidad <- ggplot() +
  geom_sf(data = peru_dep_simple, fill = "gray95", color = "gray60", size = 0.3) +
  geom_sf(data = datos_sf_intensidad, 
          aes(color = log10(intensidad_prod)), size = 1.2, alpha = 0.7) +
  scale_color_gradient2(
    low = "#2C3E50", mid = "#F39C12", high = "#E74C3C",
    midpoint = median(log10(datos_sf_intensidad$intensidad_prod), na.rm = TRUE),
    name = "Intensidad\n(log10 kg/ha)"
  ) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    title = "Intensidad de Producción Agrícola",
    subtitle = "Producción por hectárea cosechada",
    caption = "Fuente: INEI"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "right"
  )

print(mapa_intensidad)


# ============================================================================
# MAPA 3: SISTEMAS DE RIEGO
# ============================================================================

datos_sf_riego <- datos_sf %>%
  mutate(
    categoria_riego = case_when(
      prop_riego == 0 ~ "Sin riego",
      prop_riego > 0 & prop_riego < 0.5 ~ "Riego parcial",
      prop_riego >= 0.5 ~ "Con riego"
    ),
    categoria_riego = factor(categoria_riego, 
                             levels = c("Sin riego", "Riego parcial", "Con riego"))
  )

mapa_riego <- ggplot() +
  geom_sf(data = peru_dep_simple, fill = "white", color = "gray40", size = 0.4) +
  geom_sf(data = datos_sf_riego, aes(color = categoria_riego), size = 1, alpha = 0.6) +
  scale_color_manual(
    values = c("#8B4513", "#F39C12", "#27AE60"),
    name = "Sistema de Riego"
  ) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    title = "Distribución de Sistemas de Riego",
    subtitle = "Proporción de cultivos con riego",
    caption = "Fuente: INEI"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    legend.position = "right"
  )

print(mapa_riego)


# Histogramas
par(mfrow = c(2, 2))
hist(log10(datos_segmento$prod_total_kg + 1), 
     main = "Log10(Producción Total)", xlab = "Log10(kg)")
hist(datos_segmento$sup_total_ha, 
     main = "Superficie Total", xlab = "Hectáreas")
hist(log10(datos_segmento$intensidad_prod + 1), 
     main = "Log10(Intensidad de Producción)", xlab = "Log10(kg/ha)")
hist(datos_segmento$n_cultivos, 
     main = "Número de Cultivos", xlab = "N° cultivos")
par(mfrow = c(1, 1))

cat("\n=== PREPARACIÓN DE DATOS COMPLETADA ===\n")



######_-------------------------------------

# ============================================================================
# ANÁLISIS GEOESTADÍSTICO AVANZADO
# Variogramas, Kriging y Simulación Gaussiana
# ============================================================================

library(gstat)
library(sp)
library(sf)
library(tidyverse)
library(viridis)
library(automap)  # Para ajuste automático de variogramas

# Si no tienes automap:
# install.packages("automap")

cat("\n═══════════════════════════════════════════════════════\n")
cat("PARTE 2: ANÁLISIS GEOESTADÍSTICO\n")
cat("═══════════════════════════════════════════════════════\n\n")

# ============================================================================
# 1. PREPARACIÓN DE DATOS PARA GEOESTADÍSTICA
# ============================================================================

cat("1. Preparando datos para análisis geoestadístico...\n")

# Convertir datos_sf a objeto sp (necesario para gstat)
datos_sp <- as(datos_sf, "Spatial")

# Filtrar datos válidos y crear variables transformadas
datos_geo <- datos_segmento %>%
  filter(!is.na(intensidad_prod),
         !is.na(LONGITUD),
         !is.na(LATITUD),
         intensidad_prod > 0,
         intensidad_prod < quantile(intensidad_prod, 0.99, na.rm = TRUE)) %>%
  mutate(
    # Transformación Normal Score para intensidad de producción
    log_intensidad = log10(intensidad_prod + 1),
    
    # Transformación para superficie
    log_superficie = log10(sup_total_ha + 1),
    
    # Transformación para producción total
    log_produccion = log10(prod_total_kg + 1)
  )

# Convertir a objeto espacial
coordinates(datos_geo) <- ~LONGITUD+LATITUD
proj4string(datos_geo) <- CRS("+proj=longlat +datum=WGS84")

cat("   ✓ Datos preparados:", nrow(datos_geo), "observaciones\n\n")

# ============================================================================
# 2. ANÁLISIS VARIOGRÁFICO - INTENSIDAD DE PRODUCCIÓN
# ============================================================================

cat("2. Calculando variogramas empíricos...\n")

# 2.1 Variograma empírico omnidireccional
vario_omni <- variogram(log_intensidad ~ 1, 
                        data = datos_geo,
                        cutoff = 5,      # Distancia máxima en grados (~500 km)
                        width = 0.25)     # Ancho de cada lag

# 2.2 Variogramas direccionales (para detectar anisotropía)
vario_dir <- variogram(log_intensidad ~ 1,
                       data = datos_geo,
                       cutoff = 5,
                       width = 0.25,
                       alpha = c(0, 45, 90, 135))  # 4 direcciones

cat("   ✓ Variogramas calculados\n\n")

# ============================================================================
# 3. AJUSTE DE MODELOS VARIOGRÁFICOS
# ============================================================================

cat("3. Ajustando modelos teóricos al variograma...\n")

# 3.1 Ajuste automático
vario_fit_auto <- autofitVariogram(log_intensidad ~ 1, 
                                   datos_geo,
                                   model = c("Sph", "Exp", "Gau"))

# 3.2 Ajuste manual con modelo esférico
vario_model <- vgm(psill = 0.5,      # Meseta parcial
                   model = "Sph",    # Modelo esférico
                   range = 2,        # Rango
                   nugget = 0.1)     # Efecto pepita

vario_fit <- fit.variogram(vario_omni, vario_model)

cat("   ✓ Modelo ajustado:\n")
print(vario_fit)
cat("\n")

# ============================================================================
# 4. VISUALIZACIÓN DE VARIOGRAMAS
# ============================================================================

cat("4. Generando gráficos de variogramas...\n")

# 4.1 Variograma omnidireccional con modelo ajustado
p_vario_omni <- ggplot() +
  # Puntos empíricos
  geom_point(data = vario_omni, 
             aes(x = dist, y = gamma, size = np),
             color = "#3498DB", alpha = 0.6) +
  
  # Modelo ajustado
  geom_line(data = variogramLine(vario_fit, maxdist = max(vario_omni$dist)),
            aes(x = dist, y = gamma),
            color = "#E74C3C", size = 1.2) +
  
  # Línea horizontal en la meseta
  geom_hline(yintercept = sum(vario_fit$psill), 
             linetype = "dashed", color = "gray40") +
  
  # Línea vertical en el rango
  geom_vline(xintercept = vario_fit$range[2], 
             linetype = "dashed", color = "gray40") +
  
  labs(
    title = "Variograma Empírico y Modelo Teórico",
    subtitle = "Intensidad de Producción (log-transformada)",
    x = "Distancia (grados decimales)",
    y = "Semivarianza",
    size = "Pares de\npuntos",
    caption = paste0("Modelo: ", vario_fit$model[2], 
                     " | Rango: ", round(vario_fit$range[2], 2),
                     " | Meseta: ", round(sum(vario_fit$psill), 3),
                     " | Nugget: ", round(vario_fit$psill[1], 3))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray30"),
    plot.caption = element_text(size = 9, hjust = 0)
  )

print(p_vario_omni)
ggsave("variograma_omnidireccional.png", width = 10, height = 7, dpi = 300, bg = "white")

# 4.2 Variogramas direccionales (mapa de anisotropía)
p_vario_dir <- ggplot(vario_dir, aes(x = dist, y = gamma, color = factor(dir.hor))) +
  geom_point(aes(size = np), alpha = 0.6) +
  geom_line() +
  scale_color_viridis_d(
    name = "Dirección",
    labels = c("0° (N-S)", "45° (NE-SO)", "90° (E-O)", "135° (NO-SE)")
  ) +
  labs(
    title = "Variogramas Direccionales",
    subtitle = "Detección de Anisotropía Espacial",
    x = "Distancia (grados decimales)",
    y = "Semivarianza",
    size = "Pares"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p_vario_dir)
ggsave("variogramas_direccionales.png", width = 10, height = 7, dpi = 300, bg = "white")

cat("   ✓ Gráficos guardados\n\n")
#############--------------------
###--------------------------------------------




# ============================================================================
# SIMULACIÓN MULTIGAUSSIANA SECUENCIAL (SGS)
# Variables Continuas: Intensidad de Producción
# ============================================================================

library(gstat)
library(sp)
library(sf)
library(tidyverse)
library(viridis)
library(patchwork)

cat("\n═══════════════════════════════════════════════════════\n")
cat("PARTE 3: SIMULACIÓN MULTIGAUSSIANA\n")
cat("═══════════════════════════════════════════════════════\n\n")

# ============================================================================
# 1. CONFIGURACIÓN DE LA SIMULACIÓN
# ============================================================================

cat("1. Configurando parámetros de simulación...\n")

# Número de realizaciones (escenarios equiprobables)
n_simulaciones <- 100  # Puedes ajustar según capacidad: 50, 100, 200

# Grilla de predicción
lon_range <- c(-81.5, -68.5)
lat_range <- c(-18.5, -0.5)
grid_spacing <- 0.15  # ~15 km de resolución

cat("   Parámetros:\n")
cat("   • Realizaciones:", n_simulaciones, "\n")
cat("   • Resolución grilla:", grid_spacing, "grados (~", round(grid_spacing * 111, 0), "km)\n\n")

# ============================================================================
# 2. CREAR GRILLA DE SIMULACIÓN
# ============================================================================

cat("2. Creando grilla de simulación...\n")

grid_sim <- expand.grid(
  LONGITUD = seq(lon_range[1], lon_range[2], by = grid_spacing),
  LATITUD = seq(lat_range[1], lat_range[2], by = grid_spacing)
)

coordinates(grid_sim) <- ~LONGITUD+LATITUD
gridded(grid_sim) <- TRUE
proj4string(grid_sim) <- CRS("+proj=longlat +datum=WGS84")

cat("   ✓ Grilla creada:", nrow(grid_sim), "celdas\n\n")

# ============================================================================
# 3. SIMULACIÓN GAUSSIANA SECUENCIAL
# ============================================================================

cat("3. Ejecutando simulación multigaussiana...\n")
cat("   (Esto puede tomar 5-15 minutos dependiendo de tu PC)\n\n")

set.seed(12345)  # Reproducibilidad

# Ejecutar simulación con el modelo variográfico ajustado
sim_result <- krige(
  formula = log_intensidad ~ 1,
  locations = datos_geo,
  newdata = grid_sim,
  model = vario_fit,
  nmax = 30,           # Máximo de vecinos a usar
  nsim = n_simulaciones,  # Número de realizaciones
  debug.level = 0
)

cat("   ✓ Simulación completada!\n\n")

# ============================================================================
# 4. PROCESAR RESULTADOS
# ============================================================================

cat("4. Procesando resultados de simulación...\n")

# Convertir a data frame
sim_df <- as.data.frame(sim_result)

# Calcular estadísticas por celda
sim_stats <- sim_df %>%
  mutate(
    # Media de todas las realizaciones
    media = rowMeans(select(., starts_with("sim")), na.rm = TRUE),
    
    # Desviación estándar
    sd = apply(select(., starts_with("sim")), 1, sd, na.rm = TRUE),
    
    # Percentiles
    p10 = apply(select(., starts_with("sim")), 1, 
                function(x) quantile(x, 0.10, na.rm = TRUE)),
    p50 = apply(select(., starts_with("sim")), 1, 
                function(x) quantile(x, 0.50, na.rm = TRUE)),
    p90 = apply(select(., starts_with("sim")), 1, 
                function(x) quantile(x, 0.90, na.rm = TRUE)),
    
    # Coeficiente de variación (incertidumbre relativa)
    cv = sd / abs(media),
    
    # Back-transform a escala original (kg/ha)
    intensidad_media = 10^media - 1,
    intensidad_p10 = 10^p10 - 1,
    intensidad_p90 = 10^p90 - 1
  ) %>%
  select(LONGITUD, LATITUD, media, sd, cv, p10, p50, p90, 
         intensidad_media, intensidad_p10, intensidad_p90)

cat("   ✓ Estadísticas calculadas\n\n")

#########+++++++++++++++++++++++++++++++++++++++++++++++
#########+


# ============================================================================
# 7. ANÁLISIS DE PROBABILIDAD
# ============================================================================

cat("7. Calculando probabilidades de excedencia...\n")

# Definir umbrales de productividad
umbral_bajo <- 1000
umbral_medio <- 5000
umbral_alto <- 10000

# Calcular probabilidades
prob_excedencia <- sim_df %>%
  mutate(
    prob_mayor_1000 = rowMeans(
      select(., starts_with("sim")) %>%
        mutate(across(everything(), ~(10^. - 1) > umbral_bajo))
    ),
    prob_mayor_5000 = rowMeans(
      select(., starts_with("sim")) %>%
        mutate(across(everything(), ~(10^. - 1) > umbral_medio))
    ),
    prob_mayor_10000 = rowMeans(
      select(., starts_with("sim")) %>%
        mutate(across(everything(), ~(10^. - 1) > umbral_alto))
    )
  ) %>%
  select(LONGITUD, LATITUD, starts_with("prob_"))

# Mapa de probabilidad de alta productividad
mapa_prob <- ggplot() +
  geom_tile(data = prob_excedencia, 
            aes(x = LONGITUD, y = LATITUD, fill = prob_mayor_10000)) +
  geom_sf(data = peru_dep_simple, fill = NA, color = "gray30", size = 0.3) +
  scale_fill_viridis(
    option = "inferno",
    name = "Probabilidad",
    labels = scales::percent,
    na.value = "white"
  ) +
  coord_sf(xlim = lon_range, ylim = lat_range) +
  labs(
    title = "Probabilidad de Alta Productividad",
    subtitle = "P(Intensidad > 10,000 kg/ha)",
    x = "Longitud", y = "Latitud",
    caption = paste0("Basado en ", n_simulaciones, " realizaciones de SGS")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(mapa_prob)

##################
############-------


# ============================================================================
# SOLUCIÓN COMPLETA: ENMASCARAR Y RECALCULAR
# ============================================================================

cat("Aplicando máscara territorial de Perú...\n")

# 1. Convertir grilla simulada a sf
sim_sf <- st_as_sf(sim_df, coords = c("LONGITUD", "LATITUD"), crs = 4326)

# 2. Identificar qué celdas están DENTRO de Perú
celdas_en_peru <- st_intersects(sim_sf, peru_pais_simple, sparse = FALSE)

# 3. Filtrar datos ANTES de calcular estadísticas
sim_df_peru <- sim_df %>%
  mutate(dentro_peru = as.vector(celdas_en_peru)) %>%
  filter(dentro_peru == TRUE) %>%
  select(-dentro_peru)

# 4. RECALCULAR estadísticas solo para celdas en Perú
sim_stats_peru <- sim_df_peru %>%
  mutate(
    # Media
    media = rowMeans(select(., starts_with("sim")), na.rm = TRUE),
    
    # Desviación estándar
    sd = apply(select(., starts_with("sim")), 1, sd, na.rm = TRUE),
    
    # Percentiles
    p10 = apply(select(., starts_with("sim")), 1, 
                function(x) quantile(x, 0.10, na.rm = TRUE)),
    p50 = apply(select(., starts_with("sim")), 1, 
                function(x) quantile(x, 0.50, na.rm = TRUE)),
    p90 = apply(select(., starts_with("sim")), 1, 
                function(x) quantile(x, 0.90, na.rm = TRUE)),
    
    # Coeficiente de variación
    cv = sd / abs(media),
    
    # Back-transform a escala original
    intensidad_media = 10^media - 1,
    intensidad_p10 = 10^p10 - 1,
    intensidad_p90 = 10^p90 - 1,
    
    # IMPORTANTE: Calcular rango aquí
    rango_p90_p10 = intensidad_p90 - intensidad_p10
  ) %>%
  select(LONGITUD, LATITUD, media, sd, cv, p10, p50, p90, 
         intensidad_media, intensidad_p10, intensidad_p90, rango_p90_p10)

cat("   ✓ Celdas totales:", nrow(sim_df), "\n")
cat("   ✓ Celdas dentro de Perú:", nrow(sim_stats_peru), "\n\n")

# ============================================================================
# 5. VISUALIZACIONES CON MÁSCARA APLICADA
# ============================================================================

cat("5. Generando mapas SOLO dentro de Perú...\n")

# 5.1 Mapa de MEDIA
mapa_media <- ggplot() +
  geom_sf(data = peru_dep_simple, fill = "gray95", color = "gray40", size = 0.4) +
  geom_tile(data = sim_stats_peru, 
            aes(x = LONGITUD, y = LATITUD, fill = intensidad_media)) +
  geom_sf(data = peru_dep_simple, fill = NA, color = "gray20", size = 0.3) +
  scale_fill_gradientn(
    colors = c("#440154", "#31688E", "#35B779", "#FDE724"),
    name = "Intensidad\nMedia\n(kg/ha)",
    trans = "log10",
    labels = scales::comma,
    na.value = NA
  ) +
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  labs(
    title = "Simulación Multigaussiana: Intensidad de Producción",
    subtitle = paste0("Media de ", n_simulaciones, " realizaciones - Territorio peruano"),
    x = "Longitud", y = "Latitud",
    caption = "Fuente: INEI | Método: Sequential Gaussian Simulation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
    legend.position = "right"
  )

print(mapa_media)
ggsave("SGS_mapa_media_PERU.png", width = 12, height = 10, dpi = 300, bg = "white")

# 5.2 Mapa de INCERTIDUMBRE
mapa_cv <- ggplot() +
  geom_sf(data = peru_dep_simple, fill = "gray95", color = "gray40", size = 0.4) +
  geom_tile(data = sim_stats_peru %>% filter(cv < 1.5), 
            aes(x = LONGITUD, y = LATITUD, fill = cv)) +
  geom_sf(data = peru_dep_simple, fill = NA, color = "gray20", size = 0.3) +
  scale_fill_gradientn(
    colors = c("#FDE724", "#35B779", "#31688E", "#440154"),
    name = "Coeficiente\nde Variación",
    na.value = NA
  ) +
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  labs(
    title = "Mapa de Incertidumbre",
    subtitle = "Coeficiente de variación - Mayor valor = mayor incertidumbre",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right"
  )

print(mapa_cv)
ggsave("SGS_mapa_incertidumbre_PERU.png", width = 12, height = 10, dpi = 300, bg = "white")

# 5.3 Mapa de RANGO (AHORA SÍ FUNCIONARÁ)
mapa_rango <- ggplot() +
  geom_sf(data = peru_dep_simple, fill = "gray95", color = "gray40", size = 0.4) +
  geom_tile(data = sim_stats_peru, 
            aes(x = LONGITUD, y = LATITUD, fill = rango_p90_p10)) +
  geom_sf(data = peru_dep_simple, fill = NA, color = "gray20", size = 0.3) +
  scale_fill_gradientn(
    colors = c("#0D0887", "#7E03A8", "#CC4678", "#F89441", "#F0F921"),
    name = "Rango\nP90-P10\n(kg/ha)",
    trans = "log10",
    labels = scales::comma,
    na.value = NA
  ) +
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  labs(
    title = "Amplitud del Intervalo de Confianza (80%)",
    subtitle = "Diferencia entre percentil 90 y percentil 10",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right"
  )

print(mapa_rango)
ggsave("SGS_mapa_rango_PERU.png", width = 12, height = 10, dpi = 300, bg = "white")

cat("   ✓ 3 mapas guardados\n\n")

# ============================================================================
# 6. PANEL DE REALIZACIONES
# ============================================================================

cat("6. Creando panel de realizaciones...\n")

set.seed(456)
realizaciones_muestra <- sample(1:n_simulaciones, 6)

mapas_realizaciones <- list()

for(i in 1:6) {
  sim_id <- realizaciones_muestra[i]
  col_name <- paste0("sim", sim_id)
  
  datos_plot <- sim_df_peru %>%
    select(LONGITUD, LATITUD, valor = all_of(col_name)) %>%
    filter(!is.na(valor)) %>%
    mutate(intensidad = pmax(1, 10^valor - 1))
  
  p <- ggplot() +
    geom_sf(data = peru_dep_simple, fill = "gray95", color = "gray40", size = 0.2) +
    geom_tile(data = datos_plot, 
              aes(x = LONGITUD, y = LATITUD, fill = intensidad)) +
    geom_sf(data = peru_dep_simple, fill = NA, color = "gray20", size = 0.15) +
    scale_fill_gradientn(
      colors = c("#440154", "#31688E", "#35B779", "#FDE724"),
      trans = "log10",
      name = "kg/ha",
      guide = guide_colorbar(barwidth = 0.6, barheight = 5)
    ) +
    coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
    labs(title = paste("Realización", sim_id)) +
    theme_void() +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7)
    )
  
  mapas_realizaciones[[i]] <- p
}

panel_realizaciones <- wrap_plots(mapas_realizaciones, ncol = 3) +
  plot_annotation(
    title = "Realizaciones Equiprobables - Territorio Peruano",
    subtitle = "Cada mapa representa un escenario posible con igual probabilidad",
    caption = paste0("Sequential Gaussian Simulation | n=", n_simulaciones, " realizaciones"),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      plot.caption = element_text(size = 9, hjust = 0.5)
    )
  )

print(panel_realizaciones)
ggsave("SGS_panel_realizaciones_PERU.png", width = 16, height = 12, dpi = 300, bg = "white")

cat("   ✓ Panel guardado\n\n")

# ============================================================================
# 7. GUARDAR RESULTADOS
# ============================================================================

cat("7. Guardando resultados finales...\n")

write_csv(sim_stats_peru, "SGS_estadisticas_PERU.csv")
save(sim_stats_peru, sim_df_peru, vario_fit, n_simulaciones,
     file = "resultados_SGS_PERU.RData")

cat("   ✓ Archivos guardados\n\n")

cat("═══════════════════════════════════════════════════════\n")
cat("SIMULACIÓN MULTIGAUSSIANA COMPLETADA ✓\n")
cat("Solo territorio peruano | ", nrow(sim_stats_peru), " celdas\n")
cat("═══════════════════════════════════════════════════════\n")




# ============================================================================
# 8. GUARDAR RESULTADOS
# ============================================================================

cat("8. Guardando resultados...\n")

# Guardar estadísticas
write_csv(sim_stats, "SGS_estadisticas_por_celda.csv")

# Guardar workspace
save(sim_result, sim_stats, prob_excedencia, vario_fit,
     file = "resultados_simulacion_multigaussiana.RData")

cat("   ✓ Archivos guardados\n\n")

# ============================================================================
# 9. RESUMEN
# ============================================================================

cat("\n═══════════════════════════════════════════════════════\n")
cat("RESUMEN DE SIMULACIÓN MULTIGAUSSIANA\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("CONFIGURACIÓN:\n")
cat("  • Realizaciones:", n_simulaciones, "\n")
cat("  • Celdas simuladas:", nrow(sim_stats), "\n")
cat("  • Modelo variográfico: Esférico\n")
cat("  • Rango:", round(vario_fit$range[2], 2), "grados (~", 
    round(vario_fit$range[2] * 111, 0), "km)\n\n")

cat("ESTADÍSTICAS GLOBALES:\n")
cat("  • Intensidad media:", round(mean(sim_stats$intensidad_media, na.rm = TRUE), 0), "kg/ha\n")
cat("  • CV promedio:", round(mean(sim_stats$cv, na.rm = TRUE), 2), "\n")
cat("  • P10 promedio:", round(mean(sim_stats$intensidad_p10, na.rm = TRUE), 0), "kg/ha\n")
cat("  • P90 promedio:", round(mean(sim_stats$intensidad_p90, na.rm = TRUE), 0), "kg/ha\n\n")



cat("═══════════════════════════════════════════════════════\n")
cat("SIMULACIÓN MULTIGAUSSIANA COMPLETADA ✓\n")
cat("═══════════════════════════════════════════════════════\n")

########################
#----------------------------------------------------------------------------


# ============================================================================
# SIMULACIÓN PLURIGAUSSIANA
# Variables Categóricas: Tipo de Riego y Fuente de Agua
# ============================================================================

library(gstat)
library(sp)
library(sf)
library(tidyverse)
library(viridis)

cat("\n═══════════════════════════════════════════════════════\n")
cat("PARTE 4: SIMULACIÓN PLURIGAUSSIANA\n")
cat("═══════════════════════════════════════════════════════\n\n")

# ============================================================================
# 1. PREPARAR DATOS CATEGÓRICOS
# ============================================================================

cat("1. Preparando datos de irrigación...\n")

# Crear variable categórica de tipo de riego
datos_riego <- datos_segmento %>%
  filter(!is.na(LONGITUD), !is.na(LATITUD)) %>%
  mutate(
    # Categoría de riego (simplificada a 3 clases)
    tipo_riego = case_when(
      prop_riego == 0 ~ 1,                    # Sin riego
      prop_riego > 0 & prop_riego < 0.5 ~ 2,  # Riego parcial
      prop_riego >= 0.5 ~ 3                    # Con riego
    ),
    tipo_riego_nombre = case_when(
      tipo_riego == 1 ~ "Sin riego",
      tipo_riego == 2 ~ "Riego parcial",
      tipo_riego == 3 ~ "Con riego"
    )
  ) %>%
  filter(!is.na(tipo_riego))

# Convertir a objeto espacial
coordinates(datos_riego) <- ~LONGITUD+LATITUD
proj4string(datos_riego) <- CRS("+proj=longlat +datum=WGS84")

cat("   ✓ Datos preparados:", nrow(datos_riego), "observaciones\n")
cat("   ✓ Distribución de categorías:\n")
print(table(datos_riego$tipo_riego_nombre))
cat("\n")

# ============================================================================
# 2. ANÁLISIS DE AUTOCORRELACIÓN CATEGÓRICA
# ============================================================================

cat("2. Analizando autocorrelación espacial de categorías...\n")

# Crear vecindad
coords_riego <- coordinates(datos_riego)
vecinos_riego <- knearneigh(coords_riego, k = 8)
vecinos_nb <- knn2nb(vecinos_riego)

# Join-count statistics (equivalente a Moran para categorías)
join_test <- joincount.test(as.factor(datos_riego$tipo_riego),
                            nb2listw(vecinos_nb, style = "B"))

cat("   Resultados Join-Count:\n")
for(i in 1:length(join_test)) {
  cat("     ", names(join_test)[i], "\n")
  cat("       Observado:", join_test[[i]]$estimate[1], "\n")
  cat("       Esperado:", join_test[[i]]$estimate[2], "\n")
  cat("       p-valor:", format.pval(join_test[[i]]$p.value), "\n\n")
}

# ============================================================================
# 3. MÉTODO DE INDICADORES (SIMPLIFICADO)
# ============================================================================

cat("3. Aplicando método de indicadores para simulación plurigaussiana...\n\n")

# Crear indicadores binarios para cada categoría
datos_indicadores <- as.data.frame(datos_riego) %>%
  mutate(
    ind_sin_riego = as.numeric(tipo_riego == 1),
    ind_parcial = as.numeric(tipo_riego == 2),
    ind_con_riego = as.numeric(tipo_riego == 3)
  )

# Convertir a espacial
coordinates(datos_indicadores) <- ~LONGITUD+LATITUD
proj4string(datos_indicadores) <- CRS("+proj=longlat +datum=WGS84")

# ============================================================================
# 4. VARIOGRAMAS PARA CADA INDICADOR
# ============================================================================

cat("4. Calculando variogramas de indicadores...\n")

# Variograma para "Sin riego"
vario_ind1 <- variogram(ind_sin_riego ~ 1, datos_indicadores, 
                        cutoff = 5, width = 0.25)
fit_ind1 <- fit.variogram(vario_ind1, 
                          vgm(psill = 0.2, model = "Sph", range = 2, nugget = 0.05))

# Variograma para "Con riego"
vario_ind3 <- variogram(ind_con_riego ~ 1, datos_indicadores, 
                        cutoff = 5, width = 0.25)
fit_ind3 <- fit.variogram(vario_ind3, 
                          vgm(psill = 0.2, model = "Sph", range = 2, nugget = 0.05))

cat("   ✓ Variogramas ajustados para cada indicador\n\n")

# Visualizar variogramas de indicadores
p_vario_ind <- ggplot() +
  geom_point(data = vario_ind1, aes(x = dist, y = gamma, color = "Sin riego"), size = 2) +
  geom_line(data = variogramLine(fit_ind1, maxdist = 5), 
            aes(x = dist, y = gamma, color = "Sin riego"), size = 1) +
  geom_point(data = vario_ind3, aes(x = dist, y = gamma, color = "Con riego"), size = 2) +
  geom_line(data = variogramLine(fit_ind3, maxdist = 5), 
            aes(x = dist, y = gamma, color = "Con riego"), size = 1) +
  scale_color_manual(values = c("Sin riego" = "#E74C3C", "Con riego" = "#27AE60"),
                     name = "Categoría") +
  labs(
    title = "Variogramas de Indicadores - Tipo de Riego",
    x = "Distancia (grados)",
    y = "Semivarianza"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p_vario_ind)
ggsave("PLURI_variogramas_indicadores.png", width = 10, height = 7, dpi = 300, bg = "white")

# ============================================================================
# 5. SIMULACIÓN DE INDICADORES
# ============================================================================

cat("5. Ejecutando simulación de indicadores...\n")
cat("   (Generando 50 realizaciones por categoría)\n\n")

n_sim_pluri <- 50

# Usar la misma grilla que en multigaussiana, pero filtrada para Perú
grid_riego <- expand.grid(
  LONGITUD = seq(lon_range[1], lon_range[2], by = grid_spacing),
  LATITUD = seq(lat_range[1], lat_range[2], by = grid_spacing)
)

coordinates(grid_riego) <- ~LONGITUD+LATITUD
gridded(grid_riego) <- TRUE
proj4string(grid_riego) <- CRS("+proj=longlat +datum=WGS84")

set.seed(789)

# Simular indicador "Sin riego"
sim_ind1 <- krige(ind_sin_riego ~ 1,
                  locations = datos_indicadores,
                  newdata = grid_riego,
                  model = fit_ind1,
                  nmax = 30,
                  nsim = n_sim_pluri)

# Simular indicador "Con riego"
sim_ind3 <- krige(ind_con_riego ~ 1,
                  locations = datos_indicadores,
                  newdata = grid_riego,
                  model = fit_ind3,
                  nmax = 30,
                  nsim = n_sim_pluri)

cat("   ✓ Simulaciones de indicadores completadas\n\n")

# ============================================================================
# 6. ASIGNAR CATEGORÍAS (MÉTODO DE MÁXIMA PROBABILIDAD)
# ============================================================================

cat("6. Asignando categorías finales...\n")

# Convertir a data frames
df_ind1 <- as.data.frame(sim_ind1)
df_ind3 <- as.data.frame(sim_ind3)

# Calcular probabilidades promedio por celda
prob_categorias <- df_ind1 %>%
  select(LONGITUD, LATITUD) %>%
  mutate(
    # Probabilidad de "Sin riego" (promedio de simulaciones)
    prob_sin_riego = rowMeans(df_ind1[, grepl("sim", names(df_ind1))], na.rm = TRUE),
    
    # Probabilidad de "Con riego"
    prob_con_riego = rowMeans(df_ind3[, grepl("sim", names(df_ind3))], na.rm = TRUE),
    
    # Probabilidad de "Parcial" (complemento)
    prob_parcial = pmax(0, 1 - prob_sin_riego - prob_con_riego),
    
    # Normalizar para que sumen 1
    suma = prob_sin_riego + prob_parcial + prob_con_riego,
    prob_sin_riego = prob_sin_riego / suma,
    prob_parcial = prob_parcial / suma,
    prob_con_riego = prob_con_riego / suma,
    
    # Categoría más probable
    categoria_final = case_when(
      prob_sin_riego > prob_parcial & prob_sin_riego > prob_con_riego ~ "Sin riego",
      prob_con_riego > prob_sin_riego & prob_con_riego > prob_parcial ~ "Con riego",
      TRUE ~ "Riego parcial"
    )
  ) %>%
  select(-suma)

cat("   ✓ Categorías asignadas\n\n")

# ============================================================================
# 7. APLICAR MÁSCARA DE PERÚ
# ============================================================================

cat("7. Aplicando máscara territorial...\n")

prob_sf <- st_as_sf(prob_categorias, coords = c("LONGITUD", "LATITUD"), crs = 4326)
dentro_peru_pluri <- st_intersects(prob_sf, peru_pais_simple, sparse = FALSE)

prob_categorias_peru <- prob_categorias %>%
  mutate(dentro = as.vector(dentro_peru_pluri)) %>%
  filter(dentro == TRUE) %>%
  select(-dentro)

cat("   ✓ Celdas en Perú:", nrow(prob_categorias_peru), "\n\n")

# ============================================================================
# 8. VISUALIZACIONES
# ============================================================================

cat("8. Generando mapas de simulación plurigaussiana...\n")

# 8.1 Mapa de categoría más probable
mapa_categoria <- ggplot() +
  geom_sf(data = peru_dep_simple, fill = "gray95", color = "gray40", size = 0.4) +
  geom_tile(data = prob_categorias_peru,
            aes(x = LONGITUD, y = LATITUD, fill = categoria_final)) +
  geom_sf(data = peru_dep_simple, fill = NA, color = "gray20", size = 0.3) +
  scale_fill_manual(
    values = c("Sin riego" = "#8B4513",
               "Riego parcial" = "#F39C12",
               "Con riego" = "#27AE60"),
    name = "Tipo de Riego\n(Categoría más probable)"
  ) +
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  labs(
    title = "Simulación Plurigaussiana: Tipo de Riego Predominante",
    subtitle = "Método de indicadores con asignación por máxima probabilidad",
    x = "Longitud", y = "Latitud",
    caption = "Fuente: INEI | Método: Indicator Simulation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.position = "right"
  )

print(mapa_categoria)
ggsave("PLURI_mapa_categoria_final.png", width = 12, height = 10, dpi = 300, bg = "white")

# 8.2 Mapa de probabilidad de "Con riego"
mapa_prob_riego <- ggplot() +
  geom_sf(data = peru_dep_simple, fill = "gray95", color = "gray40", size = 0.4) +
  geom_tile(data = prob_categorias_peru,
            aes(x = LONGITUD, y = LATITUD, fill = prob_con_riego)) +
  geom_sf(data = peru_dep_simple, fill = NA, color = "gray20", size = 0.3) +
  scale_fill_viridis(
    option = "inferno",
    name = "Probabilidad",
    labels = scales::percent
  ) +
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  labs(
    title = "Probabilidad de Presencia de Riego",
    subtitle = "P(Tipo de Riego = 'Con riego')",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(mapa_prob_riego)
ggsave("PLURI_probabilidad_con_riego.png", width = 12, height = 10, dpi = 300, bg = "white")

# 8.3 Mapa de incertidumbre categórica (entropía)
prob_categorias_peru <- prob_categorias_peru %>%
  mutate(
    entropia = -(prob_sin_riego * log(prob_sin_riego + 1e-10) +
                   prob_parcial * log(prob_parcial + 1e-10) +
                   prob_con_riego * log(prob_con_riego + 1e-10))
  )

mapa_entropia <- ggplot() +
  geom_sf(data = peru_dep_simple, fill = "gray95", color = "gray40", size = 0.4) +
  geom_tile(data = prob_categorias_peru,
            aes(x = LONGITUD, y = LATITUD, fill = entropia)) +
  geom_sf(data = peru_dep_simple, fill = NA, color = "gray20", size = 0.3) +
  scale_fill_viridis(
    option = "magma",
    name = "Entropía\n(Incertidumbre)",
    direction = -1
  ) +
  coord_sf(xlim = lon_range, ylim = lat_range, expand = FALSE) +
  labs(
    title = "Incertidumbre Categórica: Entropía",
    subtitle = "Mayor entropía = Mayor incertidumbre en la categoría asignada",
    x = "Longitud", y = "Latitud"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  )

print(mapa_entropia)
ggsave("PLURI_mapa_entropia.png", width = 12, height = 10, dpi = 300, bg = "white")

cat("   ✓ 3 mapas guardados\n\n")

# ============================================================================
# 9. GUARDAR RESULTADOS
# ============================================================================

cat("9. Guardando resultados...\n")

write_csv(prob_categorias_peru, "PLURI_probabilidades_riego.csv")
save(prob_categorias_peru, fit_ind1, fit_ind3,
     file = "resultados_plurigaussiana.RData")

cat("   ✓ Archivos guardados\n\n")

# ============================================================================
# 10. RESUMEN
# ============================================================================

cat("\n═══════════════════════════════════════════════════════\n")
cat("RESUMEN DE SIMULACIÓN PLURIGAUSSIANA\n")
cat("═══════════════════════════════════════════════════════\n\n")

cat("CONFIGURACIÓN:\n")
cat("  • Realizaciones por indicador:", n_sim_pluri, "\n")
cat("  • Celdas en Perú:", nrow(prob_categorias_peru), "\n")
cat("  • Método: Simulación de indicadores\n\n")

cat("DISTRIBUCIÓN DE CATEGORÍAS:\n")
print(table(prob_categorias_peru$categoria_final))
cat("\n")


cat("═══════════════════════════════════════════════════════\n")
cat("SIMULACIÓN PLURIGAUSSIANA COMPLETADA ✓\n")
cat("═══════════════════════════════════════════════════════\n")
