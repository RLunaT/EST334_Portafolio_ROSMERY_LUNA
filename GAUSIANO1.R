# =============================================================================
# ANÁLISIS GAUSSIANO FUNCIONAL - DATOS AGRÍCOLAS
# =============================================================================

# Librerías necesarias
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(MASS)
library(corrplot)
library(nortest)
library(mixtools)
library(viridis)
base <- read_csv("C:/Users/ROSMERY/Downloads/base.csv") 

# =============================================================================
# 1. PREPARAR DATOS
# =============================================================================

datos <- base %>%
  filter(!is.na(P217_SUP_ha) & !is.na(P219_CANT_1)) %>%
  dplyr::select(
    superficie_ha = P217_SUP_ha,
    produccion_kg = P219_CANT_1,
    precio_kg = P220_1_PRE_KG,
    valor_total = P220_1_VAL,
    producto = P204_NOM
  ) %>%
  filter(superficie_ha > 0 & superficie_ha < 1000,
         produccion_kg > 0) %>%
  mutate(
    rendimiento = produccion_kg / superficie_ha,
    log_superficie = log(superficie_ha),
    log_produccion = log(produccion_kg),
    log_rendimiento = log(rendimiento)
  ) %>%
  filter(complete.cases(.))

cat("✓ Datos preparados:", nrow(datos), "observaciones\n")

# =============================================================================
# 2. ANÁLISIS DE NORMALIDAD
# =============================================================================

# Variables a probar
vars_test <- c("log_superficie", "log_produccion", "log_rendimiento")

# Tests simples de normalidad
resultados <- data.frame()
for(var in vars_test) {
  x <- datos[[var]]
  ad_p <- ad.test(x)$p.value
  lillie_p <- lillie.test(x)$p.value
  skew <- e1071::skewness(x)
  
  resultados <- rbind(resultados, data.frame(
    variable = var,
    anderson_p = round(ad_p, 6),
    lilliefors_p = round(lillie_p, 6),
    skewness = round(skew, 4)
  ))
}

print("=== TESTS DE NORMALIDAD ===")
print(resultados)

# Seleccionar mejores variables
var1 <- "log_superficie"
var2 <- "log_rendimiento"
cat("\n✓ Variables seleccionadas:", var1, "y", var2, "\n")

# =============================================================================
# 3. PARÁMETROS DEL CAMPO GAUSSIANO
# =============================================================================

# Matriz de datos bivariados
X <- as.matrix(datos[, c(var1, var2)])
n_obs <- nrow(X)

# Parámetros empíricos
mu_emp <- colMeans(X)
sigma_emp <- cov(X)
correlacion <- cor(datos[[var1]], datos[[var2]])

cat("\n=== PARÁMETROS DEL CAMPO GAUSSIANO ===\n")
cat("Media", var1, ":", round(mu_emp[1], 4), "\n")
cat("Media", var2, ":", round(mu_emp[2], 4), "\n")
cat("Desviación", var1, ":", round(sqrt(sigma_emp[1,1]), 4), "\n")
cat("Desviación", var2, ":", round(sqrt(sigma_emp[2,2]), 4), "\n")
cat("Correlación:", round(correlacion, 4), "\n")

# =============================================================================
# 4. SIMULACIONES DEL CAMPO GAUSSIANO
# =============================================================================

set.seed(123)

# Simulación principal
simulaciones <- MASS::mvrnorm(n = 2000, mu = mu_emp, Sigma = sigma_emp)
sim_df <- data.frame(
  log_superficie_sim = simulaciones[,1],
  log_rendimiento_sim = simulaciones[,2]
)

# Escenarios adicionales
sim_conservador <- MASS::mvrnorm(500, mu = mu_emp * 0.9, Sigma = sigma_emp * 0.8)
sim_optimista <- MASS::mvrnorm(500, mu = mu_emp * 1.1, Sigma = sigma_emp * 1.2)

cat("\n✓ Simulaciones generadas:", nrow(sim_df), "principales\n")

# =============================================================================
# 5. VISUALIZACIONES
# =============================================================================

cat("\n=== GENERANDO VISUALIZACIONES ===\n")

# 5.1 Distribuciones marginales
p1 <- ggplot(datos, aes(x = log_superficie)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.7, fill = "lightblue") +
  geom_density(color = "red", linewidth = 1) +
  stat_function(fun = dnorm, 
                args = list(mean = mu_emp[1], sd = sqrt(sigma_emp[1,1])),
                color = "blue", linewidth = 1, linetype = "dashed") +
  labs(title = "Distribución Log(Superficie)", 
       subtitle = "Observada (rojo) vs Normal teórica (azul)",
       x = "Log(Superficie)", y = "Densidad") +
  theme_minimal()

p2 <- ggplot(datos, aes(x = log_rendimiento)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, alpha = 0.7, fill = "lightgreen") +
  geom_density(color = "red", linewidth = 1) +
  stat_function(fun = dnorm, 
                args = list(mean = mu_emp[2], sd = sqrt(sigma_emp[2,2])),
                color = "blue", linewidth = 1, linetype = "dashed") +
  labs(title = "Distribución Log(Rendimiento)", 
       subtitle = "Observada (rojo) vs Normal teórica (azul)",
       x = "Log(Rendimiento)", y = "Densidad") +
  theme_minimal()

# Mostrar distribuciones
grid.arrange(p1, p2, ncol = 2)

# 5.2 Q-Q Plots
p3 <- ggplot(datos, aes(sample = log_superficie)) +
  stat_qq(alpha = 0.6) + 
  stat_qq_line(color = "red", linewidth = 1) +
  labs(title = "Q-Q Plot: Log(Superficie)") +
  theme_minimal()

p4 <- ggplot(datos, aes(sample = log_rendimiento)) +
  stat_qq(alpha = 0.6) + 
  stat_qq_line(color = "red", linewidth = 1) +
  labs(title = "Q-Q Plot: Log(Rendimiento)") +
  theme_minimal()

grid.arrange(p3, p4, ncol = 2)

# 5.3 Campo gaussiano bivariado (con manejo de errores)
tryCatch({
  # Limpiar dispositivo gráfico
  if(dev.cur() != 1) dev.off()
  
  p_campo <- ggplot() +
    geom_point(data = datos, aes(x = log_superficie, y = log_rendimiento), 
               color = "red", alpha = 0.5, size = 0.8) +
    geom_point(data = sim_df, aes(x = log_superficie_sim, y = log_rendimiento_sim), 
               color = "blue", alpha = 0.3, size = 0.6) +
    labs(title = "Campo Gaussiano Aleatorio Bivariado",
         subtitle = "Datos reales (rojo) vs Simulaciones (azul)",
         x = "Log(Superficie)", y = "Log(Rendimiento)") +
    theme_minimal()
  
  print(p_campo)
}, error = function(e) {
  cat("Error en gráfico bivariado:", e$message, "\n")
  cat("Intentando alternativa...\n")
  
  # Versión simplificada
  plot(datos$log_superficie, datos$log_rendimiento, 
       col = "red", pch = 19, cex = 0.5, alpha = 0.5,
       main = "Campo Gaussiano Bivariado",
       xlab = "Log(Superficie)", ylab = "Log(Rendimiento)")
  points(sim_df$log_superficie_sim, sim_df$log_rendimiento_sim,
         col = "blue", pch = 19, cex = 0.3)
  legend("topright", c("Datos reales", "Simulaciones"), 
         col = c("red", "blue"), pch = 19)
})

# 5.4 Contornos de densidad
p_contornos <- ggplot() +
  geom_density_2d(data = datos, aes(x = log_superficie, y = log_rendimiento), 
                  color = "red", linewidth = 1, alpha = 0.8) +
  geom_density_2d(data = sim_df, aes(x = log_superficie_sim, y = log_rendimiento_sim), 
                  color = "blue", linewidth = 1, alpha = 0.8) +
  geom_point(data = datos, aes(x = log_superficie, y = log_rendimiento),
             alpha = 0.3, size = 0.5, color = "darkred") +
  labs(title = "Contornos de Densidad del Campo Gaussiano",
       subtitle = "Datos reales (rojo) vs Simulaciones (azul)",
       x = "Log(Superficie)", y = "Log(Rendimiento)") +
  theme_minimal()

print(p_contornos)

# 5.5 Mapa de calor
x_range <- range(X[,1])
y_range <- range(X[,2])
x_grid <- seq(x_range[1], x_range[2], length.out = 50)
y_grid <- seq(y_range[1], y_range[2], length.out = 50)
grid_points <- expand.grid(x = x_grid, y = y_grid)

# Calcular densidad gaussiana (usando mixtools)
densidad <- apply(grid_points, 1, function(punto) {
  mixtools::dmvnorm(punto, mu = mu_emp, sigma = sigma_emp)
})

grid_df <- data.frame(
  x = grid_points$x,
  y = grid_points$y,
  densidad = densidad
)

p_heatmap <- ggplot(grid_df, aes(x = x, y = y, fill = densidad)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Densidad", option = "plasma") +
  geom_point(data = datos, aes(x = log_superficie, y = log_rendimiento), 
             inherit.aes = FALSE, alpha = 0.4, size = 0.5, color = "white") +
  labs(title = "Campo Gaussiano - Mapa de Densidad",
       subtitle = "Puntos blancos: observaciones",
       x = "Log(Superficie)", y = "Log(Rendimiento)") +
  theme_minimal()

print(p_heatmap)
# =============================================================================
# 6. ESCENARIOS COMPARATIVOS
# =============================================================================

# Crear dataframe combinado de escenarios
escenarios_df <- bind_rows(
  data.frame(sim_df, escenario = "Base"),
  data.frame(log_superficie_sim = sim_conservador[,1], 
             log_rendimiento_sim = sim_conservador[,2], 
             escenario = "Conservador"),
  data.frame(log_superficie_sim = sim_optimista[,1], 
             log_rendimiento_sim = sim_optimista[,2], 
             escenario = "Optimista")
)

names(escenarios_df)[1:2] <- c("log_superficie", "log_rendimiento")

p_escenarios <- ggplot(escenarios_df, aes(x = log_superficie, y = log_rendimiento, color = escenario)) +
  geom_point(alpha = 0.6, size = 0.8) +
  stat_ellipse(level = 0.95, linewidth = 1.2) +
  scale_color_manual(values = c("Base" = "#1f78b4", 
                                "Conservador" = "#e31a1c", 
                                "Optimista" = "#33a02c")) +
  labs(title = "Simulación de Escenarios del Campo Gaussiano",
       subtitle = "Diferentes parametrizaciones del modelo",
       x = "Log(Superficie)", y = "Log(Rendimiento)",
       color = "Escenario") +
  theme_minimal()

print(p_escenarios)

# =============================================================================
# 7. MATRIZ DE CORRELACIÓN
# =============================================================================

# Correlaciones entre variables originales
vars_corr <- datos %>% 
  dplyr::select(superficie_ha, produccion_kg, rendimiento, precio_kg) %>%
  dplyr::select(where(is.numeric))

if(ncol(vars_corr) >= 2) {
  matriz_corr <- cor(vars_corr, use = "complete.obs")
  corrplot(matriz_corr, method = "color", type = "upper",
           title = "Matriz de Correlación - Variables Principales", 
           mar = c(0,0,2,0))
}

# =============================================================================
# 8. VALIDACIÓN
# =============================================================================

# Test de bondad de ajuste
sim_test <- sim_df[1:nrow(datos), ]
ks_test1 <- ks.test(datos$log_superficie, sim_test$log_superficie_sim)
ks_test2 <- ks.test(datos$log_rendimiento, sim_test$log_rendimiento_sim)

cat("\n=== VALIDACIÓN DEL MODELO ===\n")
cat("Test Kolmogorov-Smirnov:\n")
cat(" • log_superficie - p-valor:", round(ks_test1$p.value, 6), "\n")
cat(" • log_rendimiento - p-valor:", round(ks_test2$p.value, 6), "\n")

# =============================================================================
# 9. RESUMEN FINAL
# =============================================================================

cat("\n" %+% strrep("=", 50) %+% "\n")
cat("RESUMEN FINAL - ANÁLISIS GAUSSIANO\n")
cat(strrep("=", 50) %+% "\n")
cat("DATOS PROCESADOS:\n")
cat(" • Observaciones analizadas:", nrow(datos), "\n")
cat(" • Variables del campo:", var1, "y", var2, "\n")

cat("\nPARÁMETROS DEL CAMPO GAUSSIANO:\n")
cat(" • Media log_superficie:", round(mu_emp[1], 4), "\n")
cat(" • Media log_rendimiento:", round(mu_emp[2], 4), "\n")
cat(" • Correlación bivariada:", round(correlacion, 4), "\n")

cat("\nSIMULACIONES:\n")
cat(" • Realizaciones principales:", nrow(sim_df), "\n")
cat(" • Escenarios generados: 3 (base, conservador, optimista)\n")

cat("\nVALIDACIÓN:\n")
cat(" • Tests KS p-valores:", round(ks_test1$p.value, 6), "y", round(ks_test2$p.value, 6), "\n")

# Guardar resultados
resultados_finales <- list(
  datos = datos,
  parametros = list(media = mu_emp, covarianza = sigma_emp, correlacion = correlacion),
  simulaciones = sim_df,
  escenarios = escenarios_df,
  validacion = list(ks1 = ks_test1$p.value, ks2 = ks_test2$p.value)
)

cat("\n=== ANÁLISIS COMPLETADO ===\n")
cat("Resultados guardados en 'resultados_finales'\n")
cat("Todas las visualizaciones han sido generadas\n")

# Operador para concatenar strings
`%+%` <- function(a, b) paste0(a, b)

