# ============================================================================
# ANÁLISIS ESPACIAL COMPLETO - 8 GRÁFICOS PROFESIONALES
# [1] Matrices de Pesos, [2] Moran, [3] Geary, [4] Hotspots
# ============================================================================

# LIBRERÍAS
library(haven)
library(tidyverse)
library(sf)
library(spdep)
library(tmap)
library(viridis)
library(RColorBrewer)
library(geodata)

# ============================================================================
# 1. CARGAR DATOS
# ============================================================================

ruta_cultivos <- "C:/Users/ROSMERY/Downloads/2024/2024 - DESCOMPRIMIDO/973-Modulo1893/01_CAP100A_04.sav"
ruta_caratula <- "C:/Users/ROSMERY/Downloads/2024/2024 - DESCOMPRIMIDO/973-Modulo1893/CARATULA.sav"

cat("═══════════════════════════════════════════════\n")
cat("   ANÁLISIS ESPACIAL - PERÚ 2024\n")
cat("═══════════════════════════════════════════════\n\n")

cat("Cargando datos...\n")
cultivos <- read_sav(ruta_cultivos)
caratula <- read_sav(ruta_caratula)

# ============================================================================
# 2. PREPARAR DATOS
# ============================================================================

datos_espaciales <- cultivos %>%
  left_join(caratula, by = c("ANIO", "CCDD", "NOMBREDD", "CCPP", "NOMBREPV", 
                             "CCDI", "NOMBREDI", "NSEGM", "ID_PROD", "UA", 
                             "RESFIN", "REGION", "ESTRATO", "FACTOR_PRODUCTOR", 
                             "CODIGO")) %>%
  filter(!is.na(LATITUD) & !is.na(LONGITUD) & !is.na(P117_SUP_ha))

datos_agregados <- datos_espaciales %>%
  group_by(NSEGM, ID_PROD, LATITUD, LONGITUD, NOMBREDD, NOMBREPV, NOMBREDI, CCDD) %>%
  summarise(
    superficie_total = sum(P117_SUP_ha, na.rm = TRUE),
    n_cultivos = n(),
    .groups = "drop"
  ) %>%
  filter(superficie_total > 0)

datos_departamento <- datos_agregados %>%
  group_by(CCDD, NOMBREDD) %>%
  summarise(
    lat = mean(LATITUD, na.rm = TRUE),
    lon = mean(LONGITUD, na.rm = TRUE),
    superficie_total = sum(superficie_total, na.rm = TRUE),
    n_productores = n(),
    .groups = "drop"
  )

departamentos_sf <- st_as_sf(datos_departamento, 
                             coords = c("lon", "lat"), 
                             crs = 4326)

cat("✓ Datos preparados:", nrow(departamentos_sf), "departamentos\n")

# ============================================================================
# 3. OBTENER MAPA DE PERÚ
# ============================================================================

cat("\nDescargando límites de Perú...\n")

peru_gadm <- gadm(country = "PER", level = 1, path = tempdir())
peru <- st_as_sf(peru_gadm)

peru <- peru %>%
  mutate(NOMBREDD_clean = toupper(NAME_1))

datos_departamento <- datos_departamento %>%
  mutate(NOMBREDD_clean = toupper(NOMBREDD))

peru_datos <- peru %>%
  left_join(datos_departamento, by = "NOMBREDD_clean")

cat("✓ Mapa cargado\n")

# ============================================================================
# 4. [1] MATRICES DE PESOS ESPACIALES
# ============================================================================

cat("\n═══════════════════════════════════════════════\n")
cat("   [1] MATRICES DE PESOS ESPACIALES\n")
cat("═══════════════════════════════════════════════\n")

coords_dept <- st_coordinates(departamentos_sf)

# Matriz K-Nearest Neighbors
k <- min(5, nrow(departamentos_sf) - 1)
vecinos_knn <- knearneigh(coords_dept, k = k)
nb_knn <- knn2nb(vecinos_knn)
pesos_knn <- nb2listw(nb_knn, style = "W")

# Matriz por Distancia
dist_threshold <- 5
nb_dist <- dnearneigh(coords_dept, 0, dist_threshold)
pesos_dist <- nb2listw(nb_dist, style = "W", zero.policy = TRUE)

cat("✓ K-Nearest Neighbors (k =", k, ")\n")
cat("✓ Distancia fija (~", round(dist_threshold*111), "km)\n")

# ============================================================================
# GRÁFICO 1: MATRIZ DE CONECTIVIDAD K-NN
# ============================================================================

cat("\n─────────────────────────────────────────\n")
cat("  GRÁFICO 1: MATRIZ DE CONECTIVIDAD K-NN\n")
cat("─────────────────────────────────────────\n")

tmap_mode("plot")

# Crear líneas de conectividad
conectividad_lines <- list()
for(i in 1:length(nb_knn)) {
  vecinos <- nb_knn[[i]]
  for(j in vecinos) {
    line <- st_linestring(rbind(coords_dept[i,], coords_dept[j,]))
    conectividad_lines <- c(conectividad_lines, list(line))
  }
}

conexiones_sf <- st_sf(geometry = st_sfc(conectividad_lines, crs = 4326))

mapa_conectividad <- tm_shape(peru) +
  tm_borders(col = "gray70", lwd = 1) +
  tm_fill(col = "white") +
  tm_shape(conexiones_sf) +
  tm_lines(col = "red", lwd = 1.5, alpha = 0.6) +
  tm_shape(departamentos_sf) +
  tm_dots(size = 0.8, col = "navy", shape = 19) +
  tm_layout(
    main.title = "[1] Matriz de Pesos Espaciales - K-Nearest Neighbors (k=5)",
    main.title.size = 1.2,
    main.title.fontface = "bold",
    frame = TRUE
  ) +
  tm_compass(position = c("left", "top"), size = 2) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_credits("Cada departamento conectado a sus 5 vecinos más cercanos",
             position = c("right", "bottom"))

print(mapa_conectividad)
cat("✓ Mostrado. Presiona Enter...\n")
readline()

# ============================================================================
# 5. [2] ÍNDICE I DE MORAN
# ============================================================================

cat("\n═══════════════════════════════════════════════\n")
cat("   [2] ÍNDICE I DE MORAN\n")
cat("═══════════════════════════════════════════════\n")

moran_knn <- moran.test(departamentos_sf$superficie_total, pesos_knn)
moran_dist <- moran.test(departamentos_sf$superficie_total, pesos_dist, zero.policy = TRUE)

cat("\nI DE MORAN (K-NN):\n")
cat("  Valor observado:", round(moran_knn$estimate[1], 4), "\n")
cat("  Valor esperado:", round(moran_knn$estimate[2], 4), "\n")
cat("  p-valor:", format.pval(moran_knn$p.value, 4), "\n")
cat("  → Interpretación:", 
    ifelse(moran_knn$p.value < 0.05, 
           "★★★ AUTOCORRELACIÓN ESPACIAL SIGNIFICATIVA ★★★",
           "NO significativo"), "\n")

# ============================================================================
# GRÁFICO 2: DIAGRAMA DE DISPERSIÓN DE MORAN
# ============================================================================

cat("\n─────────────────────────────────────────\n")
cat("  GRÁFICO 2: DIAGRAMA DE MORAN\n")
cat("─────────────────────────────────────────\n")

par(mfrow = c(1, 1), mar = c(5, 5, 5, 2))

sup_std <- scale(departamentos_sf$superficie_total)[,1]
lag_sup <- lag.listw(pesos_knn, departamentos_sf$superficie_total)
lag_sup_std <- scale(lag_sup)[,1]

cuadrante <- case_when(
  sup_std > 0 & lag_sup_std > 0 ~ "High-High",
  sup_std < 0 & lag_sup_std < 0 ~ "Low-Low",
  sup_std > 0 & lag_sup_std < 0 ~ "High-Low",
  sup_std < 0 & lag_sup_std > 0 ~ "Low-High"
)

colores_cuad <- c("High-High" = "#d73027", "Low-Low" = "#4575b4",
                  "High-Low" = "#fee090", "Low-High" = "#91bfdb")

plot(sup_std, lag_sup_std,
     pch = 21, cex = 3, lwd = 2,
     col = "black",
     bg = colores_cuad[cuadrante],
     xlab = "Superficie Total (Estandarizado)",
     ylab = "Rezago Espacial (Promedio de Vecinos)",
     main = paste0("[2] Diagrama de Dispersión de Moran\n",
                   "I = ", round(moran_knn$estimate[1], 4),
                   " | p-valor = ", format.pval(moran_knn$p.value, 4)),
     cex.lab = 1.4, cex.main = 1.3)

abline(h = 0, v = 0, col = "gray40", lty = 2, lwd = 2)
abline(lm(lag_sup_std ~ sup_std), col = "red", lwd = 4)

legend("topleft", 
       legend = c("High-High (HH)", "Low-Low (LL)", "High-Low (HL)", "Low-High (LH)"),
       pt.bg = colores_cuad,
       pch = 21,
       pt.cex = 2.5,
       cex = 1.2,
       title = "Cuadrante",
       bty = "n")

outliers_idx <- abs(sup_std) > 1.5 | abs(lag_sup_std) > 1.5
if(sum(outliers_idx) > 0) {
  text(sup_std[outliers_idx], lag_sup_std[outliers_idx],
       labels = departamentos_sf$NOMBREDD[outliers_idx],
       pos = 3, cex = 0.9, font = 2, col = "black")
}

grid(col = "gray80", lty = 1)

cat("✓ Mostrado. Presiona Enter...\n")
readline()

# ============================================================================
# 6. [3] ÍNDICE C DE GEARY
# ============================================================================

cat("\n═══════════════════════════════════════════════\n")
cat("   [3] ÍNDICE C DE GEARY\n")
cat("═══════════════════════════════════════════════\n")

geary_knn <- geary.test(departamentos_sf$superficie_total, pesos_knn)
geary_dist <- geary.test(departamentos_sf$superficie_total, pesos_dist, zero.policy = TRUE)

cat("\nC DE GEARY (K-NN):\n")
cat("  Valor observado:", round(geary_knn$estimate[1], 4), "\n")
cat("  Valor esperado:", round(geary_knn$estimate[2], 4), "\n")
cat("  p-valor:", format.pval(geary_knn$p.value, 4), "\n")
cat("  → Interpretación:", 
    ifelse(geary_knn$estimate[1] < 1, "Autocorrelación POSITIVA",
           ifelse(geary_knn$estimate[1] > 1, "Autocorrelación NEGATIVA", "Aleatoria")), "\n")
cat("\n  Nota: C < 1 = Autocorrelación positiva (similar)\n")
cat("        C ≈ 1 = Distribución aleatoria\n")
cat("        C > 1 = Autocorrelación negativa (disimilar)\n")

# ============================================================================
# GRÁFICO 3: COMPARACIÓN MORAN vs GEARY
# ============================================================================

cat("\n─────────────────────────────────────────\n")
cat("  GRÁFICO 3: COMPARACIÓN MORAN vs GEARY\n")
cat("─────────────────────────────────────────\n")

par(mfrow = c(1, 2), mar = c(6, 5, 5, 2))

# Panel A: I de Moran
barplot(c(moran_knn$estimate[1], moran_knn$estimate[2]),
        names.arg = c("Observado", "Esperado\n(aleatorio)"),
        col = c("#d73027", "#fee090"),
        ylim = c(-0.1, 0.3),
        main = "[3A] Índice I de Moran",
        ylab = "I de Moran",
        cex.main = 1.3,
        cex.lab = 1.2,
        cex.names = 1.1,
        border = "black",
        lwd = 2)
abline(h = 0, col = "black", lty = 2, lwd = 2)
text(1.5, 0.25, 
     paste0("I = ", round(moran_knn$estimate[1], 4), "\n",
            "p = ", format.pval(moran_knn$p.value, 4)),
     cex = 1.2, font = 2)
if(moran_knn$p.value < 0.05) {
  text(1.5, 0.20, "★ SIGNIFICATIVO ★", col = "red", cex = 1.1, font = 2)
}

# Panel B: C de Geary
barplot(c(geary_knn$estimate[1], geary_knn$estimate[2]),
        names.arg = c("Observado", "Esperado\n(aleatorio)"),
        col = c("#4575b4", "#fee090"),
        ylim = c(0, 1.5),
        main = "[3B] Índice C de Geary",
        ylab = "C de Geary",
        cex.main = 1.3,
        cex.lab = 1.2,
        cex.names = 1.1,
        border = "black",
        lwd = 2)
abline(h = 1, col = "red", lty = 2, lwd = 2)
text(1.5, 1.3, 
     paste0("C = ", round(geary_knn$estimate[1], 4), "\n",
            "p = ", format.pval(geary_knn$p.value, 4)),
     cex = 1.2, font = 2)
text(1.5, 0.3, "C < 1 = Similar\nC > 1 = Disimilar", cex = 0.9)

cat("✓ Mostrado. Presiona Enter...\n")
readline()

# ============================================================================
# 7. [4] ANÁLISIS LISA - HOTSPOTS
# ============================================================================

cat("\n═══════════════════════════════════════════════\n")
cat("   [4] ANÁLISIS LISA - HOTSPOTS\n")
cat("═══════════════════════════════════════════════\n")

lisa <- localmoran(departamentos_sf$superficie_total, pesos_knn)

departamentos_sf <- departamentos_sf %>%
  mutate(
    lisa_I = lisa[, 1],
    lisa_pvalue = lisa[, 5],
    superficie_std = scale(superficie_total)[,1],
    lag_superficie = lag.listw(pesos_knn, superficie_total),
    lag_superficie_std = scale(lag_superficie)[,1],
    cluster_type = case_when(
      lisa_pvalue >= 0.05 ~ "No significativo",
      superficie_std > 0 & lag_superficie_std > 0 ~ "High-High (Hotspot)",
      superficie_std < 0 & lag_superficie_std < 0 ~ "Low-Low (Coldspot)",
      superficie_std > 0 & lag_superficie_std < 0 ~ "High-Low (Outlier)",
      superficie_std < 0 & lag_superficie_std > 0 ~ "Low-High (Outlier)",
      TRUE ~ "No clasificado"
    ),
    sig_categoria = case_when(
      lisa_pvalue < 0.01 ~ "p < 0.01",
      lisa_pvalue < 0.05 ~ "p < 0.05",
      lisa_pvalue < 0.10 ~ "p < 0.10",
      TRUE ~ "No sig."
    ),
    sig_categoria = factor(sig_categoria, 
                           levels = c("p < 0.01", "p < 0.05", "p < 0.10", "No sig."))
  )

tabla_clusters <- table(departamentos_sf$cluster_type)
cat("\nDistribución de Clusters LISA:\n")
print(tabla_clusters)

# Agregar al mapa
datos_dept_for_join <- departamentos_sf %>% 
  st_drop_geometry() %>%
  mutate(NOMBREDD_clean = toupper(NOMBREDD))

peru_datos <- peru_datos %>%
  left_join(datos_dept_for_join, by = "NOMBREDD_clean", suffix = c("", ".y"))

# ============================================================================
# GRÁFICO 4: MAPA DE SUPERFICIE TOTAL
# ============================================================================

cat("\n─────────────────────────────────────────\n")
cat("  GRÁFICO 4: SUPERFICIE TOTAL\n")
cat("─────────────────────────────────────────\n")

tmap_mode("plot")

mapa1 <- tm_shape(peru_datos) +
  tm_polygons("superficie_total",
              palette = "YlOrRd",
              title = "Superficie (ha)",
              border.col = "white",
              border.lwd = 1.5,
              style = "quantile",
              n = 5) +
  tm_layout(
    main.title = "[4] Superficie Total Cultivada por Departamento",
    main.title.size = 1.2,
    main.title.fontface = "bold",
    legend.outside = TRUE,
    frame = TRUE
  ) +
  tm_compass(position = c("left", "top"), size = 2) +
  tm_scale_bar(position = c("left", "bottom"))

print(mapa1)
cat("✓ Mostrado. Presiona Enter...\n")
readline()

# ============================================================================
# GRÁFICO 5: HOTSPOTS LISA
# ============================================================================

cat("\n─────────────────────────────────────────\n")
cat("  GRÁFICO 5: HOTSPOTS LISA\n")
cat("─────────────────────────────────────────\n")

colores_lisa <- c(
  "High-High (Hotspot)" = "#d73027",
  "Low-Low (Coldspot)" = "#4575b4",
  "High-Low (Outlier)" = "#fee090",
  "Low-High (Outlier)" = "#91bfdb",
  "No significativo" = "#f0f0f0"
)

mapa2 <- tm_shape(peru_datos) +
  tm_polygons("cluster_type",
              palette = colores_lisa,
              title = "Cluster",
              border.col = "gray30",
              border.lwd = 1.5) +
  tm_layout(
    main.title = "[4] Análisis de Hotspots - LISA",
    main.title.size = 1.2,
    main.title.fontface = "bold",
    legend.outside = TRUE,
    frame = TRUE
  ) +
  tm_compass(position = c("left", "top"), size = 2) +
  tm_scale_bar(position = c("left", "bottom"))

print(mapa2)
cat("✓ Mostrado. Presiona Enter...\n")
readline()

# ============================================================================
# GRÁFICO 6: SIGNIFICANCIA ESTADÍSTICA
# ============================================================================

cat("\n─────────────────────────────────────────\n")
cat("  GRÁFICO 6: SIGNIFICANCIA\n")
cat("─────────────────────────────────────────\n")

mapa3 <- tm_shape(peru_datos) +
  tm_polygons("sig_categoria",
              palette = c("#d73027", "#fc8d59", "#fee090", "#f0f0f0"),
              title = "Significancia",
              border.col = "gray30",
              border.lwd = 1.5) +
  tm_layout(
    main.title = "[4] Significancia Estadística - Test LISA",
    main.title.size = 1.2,
    main.title.fontface = "bold",
    legend.outside = TRUE,
    frame = TRUE
  ) +
  tm_compass(position = c("left", "top"), size = 2) +
  tm_scale_bar(position = c("left", "bottom"))

print(mapa3)
cat("✓ Mostrado. Presiona Enter...\n")
readline()

# ============================================================================
# GRÁFICO 7: ÍNDICE LISA LOCAL
# ============================================================================

cat("\n─────────────────────────────────────────\n")
cat("  GRÁFICO 7: ÍNDICE LISA LOCAL\n")
cat("─────────────────────────────────────────\n")

mapa4 <- tm_shape(peru_datos) +
  tm_polygons("lisa_I",
              palette = "RdBu",
              title = "I Local",
              border.col = "gray30",
              border.lwd = 1.5,
              midpoint = 0,
              style = "cont") +
  tm_layout(
    main.title = "[2/4] Índice LISA por Departamento (Moran Local)",
    main.title.size = 1.2,
    main.title.fontface = "bold",
    legend.outside = TRUE,
    frame = TRUE
  ) +
  tm_compass(position = c("left", "top"), size = 2) +
  tm_scale_bar(position = c("left", "bottom"))

print(mapa4)
cat("✓ Mostrado. Presiona Enter...\n")
readline()

# ============================================================================
# GRÁFICO 8: CORRELOGRAMA ESPACIAL
# ============================================================================

cat("\n─────────────────────────────────────────\n")
cat("  GRÁFICO 8: CORRELOGRAMA\n")
cat("─────────────────────────────────────────\n")

par(mfrow = c(1, 1), mar = c(5, 5, 5, 2))

sp_corr <- sp.correlogram(nb_knn, 
                          departamentos_sf$superficie_total,
                          order = min(5, length(nb_knn) - 1),
                          method = "I",
                          style = "W")

orden <- 1:length(sp_corr$res[,1])
valores_I <- sp_corr$res[,1]
se <- sqrt(sp_corr$res[,3])

plot(orden, valores_I,
     type = "b",
     pch = 21,
     bg = "darkblue",
     col = "darkblue",
     cex = 2.5,
     lwd = 3,
     ylim = c(min(valores_I - 2*se), max(valores_I + 2*se)),
     xlab = "Orden de Vecindad",
     ylab = "I de Moran",
     main = "[2] Correlograma Espacial - Índice I de Moran",
     cex.lab = 1.4, cex.main = 1.3)

arrows(orden, valores_I - 1.96*se, 
       orden, valores_I + 1.96*se,
       angle = 90, code = 3, length = 0.1, lwd = 2.5, col = "gray40")

abline(h = 0, col = "red", lty = 2, lwd = 3)
abline(h = -1/(length(nb_knn)-1), col = "orange", lty = 3, lwd = 2)

grid(col = "gray80", lty = 1)

legend("topright",
       legend = c("I de Moran", "IC 95%", "Aleatorio (H₀)", "Expectativa E(I)"),
       col = c("darkblue", "gray40", "red", "orange"),
       lty = c(1, 1, 2, 3),
       lwd = c(3, 2.5, 3, 2),
       pch = c(21, NA, NA, NA),
       pt.bg = c("darkblue", NA, NA, NA),
       pt.cex = 2,
       cex = 1.1,
       bty = "n")

cat("✓ Mostrado. Presiona Enter para resumen...\n")
readline()

# ============================================================================
# RESUMEN FINAL
# ============================================================================

cat("\n")
cat("═══════════════════════════════════════════════\n")
cat("   RESUMEN EJECUTIVO - ANÁLISIS ESPACIAL\n")
cat("═══════════════════════════════════════════════\n\n")

cat("DATOS ANALIZADOS:\n")
cat("  • Productores:", format(sum(datos_departamento$n_productores), big.mark = ","), "\n")
cat("  • Departamentos:", nrow(departamentos_sf), "\n")
cat("  • Superficie total:", format(round(sum(departamentos_sf$superficie_total)), 
                                    big.mark = ","), "ha\n\n")

cat("───────────────────────────────────────────────\n")
cat("[1] MATRICES DE PESOS ESPACIALES:\n")
cat("───────────────────────────────────────────────\n")
cat("  ✓ K-Nearest Neighbors (k = 5)\n")
cat("  ✓ Cada departamento tiene 5 vecinos\n")
cat("  ✓ Gráfico 1: Matriz de Conectividad\n\n")

cat("───────────────────────────────────────────────\n")
cat("[2] ÍNDICE I DE MORAN:\n")
cat("───────────────────────────────────────────────\n")
cat("  • I observado:", round(moran_knn$estimate[1], 4), "\n")
cat("  • I esperado:", round(moran_knn$estimate[2], 4), "\n")
cat("  • p-valor:", format.pval(moran_knn$p.value, 4), "\n")
cat("  • Interpretación:", 
    ifelse(moran_knn$p.value < 0.05, 
           "\n    ★★★ AUTOCORRELACIÓN ESPACIAL POSITIVA SIGNIFICATIVA ★★★",
           "NO significativo"), "\n")
cat("  • Gráficos: 2 (Diagrama), 7 (LISA Local), 8 (Correlograma)\n\n")

cat("───────────────────────────────────────────────\n")
cat("[3] ÍNDICE C DE GEARY:\n")
cat("───────────────────────────────────────────────\n")
cat("  • C observado:", round(geary_knn$estimate[1], 4), "\n")
cat("  • C esperado:", round(geary_knn$estimate[2], 4), "\n")
cat("  • p-valor:", format.pval(geary_knn$p.value, 4), "\n")
cat("  • Interpretación:", 
    ifelse(geary_knn$estimate[1] < 1, 
           "Autocorrelación POSITIVA (valores similares próximos)",
           "Distribución ALEATORIA"), "\n")
cat("  • Gráfico 3: Comparación Moran vs Geary\n\n")

cat("───────────────────────────────────────────────\n")
cat("[4] HOTSPOTS - ANÁLISIS LISA:\n")
cat("───────────────────────────────────────────────\n")
for(tipo in names(tabla_clusters)) {
  cat("  •", tipo, ":", tabla_clusters[tipo], "\n")
}

cat("\n  DEPARTAMENTOS HOTSPOT (High-High):\n")
hotspots <- departamentos_sf %>% 
  filter(cluster_type == "High-High (Hotspot)") %>%
  arrange(desc(superficie_total))

if(nrow(hotspots) > 0) {
  for(i in 1:nrow(hotspots)) {
    cat("    ", i, ".", hotspots$NOMBREDD[i], "-",
        format(round(hotspots$superficie_total[i]), big.mark = ","), "ha\n")
  }
} else {
  cat("    No se identificaron hotspots significativos\n")
}

cat("\n  • Gráficos: 4 (Superficie), 5 (Hotspots), 6 (Significancia), 7 (LISA)\n\n")

cat("═══════════════════════════════════════════════\n")
cat("   ✓ ANÁLISIS COMPLETADO - 8 GRÁFICOS\n")
cat("═══════════════════════════════════════════════\n")

# Guardar
resultados <- departamentos_sf %>%
  st_drop_geometry() %>%
  select(NOMBREDD, superficie_total, n_productores, cluster_type, lisa_I, lisa_pvalue)

write.csv(resultados, "Resultados_Analisis_Espacial_Completo.csv", row.names = FALSE)


