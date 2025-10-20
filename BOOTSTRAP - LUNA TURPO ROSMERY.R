# Script de Bootstrap en R
# Usando el dataset mtcars para diferentes análisis bootstrap

# Cargar librerías necesarias
library(ggplot2)
library(dplyr)
library(boot)
library(gridExtra)

# Cargar datos

data(mtcars)
cat("Dataset mtcars cargado:\n")
print(head(mtcars))
cat("\nDimensiones:", dim(mtcars), "\n\n")

# =============================================================================
# 1. BOOTSTRAP PARA LA MEDIA
# =============================================================================

# Función para calcular la media (para usar con boot())
media_mpg <- function(data, indices) {
  return(mean(data[indices, "mpg"]))
}

# Realizar bootstrap para la media de mpg
set.seed(123)
n_bootstrap <- 1000

bootstrap_media <- boot(data = mtcars, 
                        statistic = media_mpg, 
                        R = n_bootstrap)

print("=== BOOTSTRAP PARA LA MEDIA DE MPG ===")
print(bootstrap_media)

# Intervalo de confianza
ic_media <- boot.ci(bootstrap_media, type = c("norm", "basic", "perc", "bca"))
print("Intervalos de confianza:")
print(ic_media)

# =============================================================================
# 2. BOOTSTRAP MANUAL PASO A PASO
# =============================================================================

# Bootstrap manual para entender el proceso
set.seed(123)
n_orig <- nrow(mtcars)
medias_bootstrap <- numeric(n_bootstrap)

for(i in 1:n_bootstrap) {
  # Muestreo con reposición
  indices_muestra <- sample(1:n_orig, n_orig, replace = TRUE)
  muestra_bootstrap <- mtcars[indices_muestra, ]
  
  # Calcular estadística de interés
  medias_bootstrap[i] <- mean(muestra_bootstrap$mpg)
}

cat("\n=== BOOTSTRAP MANUAL ===\n")
cat("Media original:", mean(mtcars$mpg), "\n")
cat("Media de medias bootstrap:", mean(medias_bootstrap), "\n")
cat("Error estándar bootstrap:", sd(medias_bootstrap), "\n")
cat("IC 95% (percentiles):", quantile(medias_bootstrap, c(0.025, 0.975)), "\n")

# =============================================================================
# 3. BOOTSTRAP PARA REGRESIÓN
# =============================================================================

# Función para coeficientes de regresión
coef_regresion <- function(data, indices) {
  modelo <- lm(mpg ~ wt + hp, data = data[indices, ])
  return(coef(modelo))
}

# Bootstrap para regresión
bootstrap_reg <- boot(data = mtcars, 
                      statistic = coef_regresion, 
                      R = n_bootstrap)

print("\n=== BOOTSTRAP PARA REGRESIÓN ===")
print(bootstrap_reg)

# Intervalos de confianza para coeficientes
ic_intercepto <- boot.ci(bootstrap_reg, type = "perc", index = 1)
ic_wt <- boot.ci(bootstrap_reg, type = "perc", index = 2)
ic_hp <- boot.ci(bootstrap_reg, type = "perc", index = 3)

cat("\nIntervalos de confianza para coeficientes:\n")
cat("Intercepto:", ic_intercepto$percent[4:5], "\n")
cat("Weight (wt):", ic_wt$percent[4:5], "\n")
cat("Horsepower (hp):", ic_hp$percent[4:5], "\n")

# =============================================================================
# 4. BOOTSTRAP PARA CORRELACIÓN
# =============================================================================

# Función para correlación
correlacion_mpg_wt <- function(data, indices) {
  return(cor(data[indices, "mpg"], data[indices, "wt"]))
}

bootstrap_cor <- boot(data = mtcars, 
                      statistic = correlacion_mpg_wt, 
                      R = n_bootstrap)

print("\n=== BOOTSTRAP PARA CORRELACIÓN MPG vs WT ===")
print(bootstrap_cor)

ic_cor <- boot.ci(bootstrap_cor, type = "perc")
cat("IC 95% para correlación:", ic_cor$percent[4:5], "\n")

# =============================================================================
# 5. DEFINIR MODELO ORIGINAL PRIMERO
# =============================================================================

# Definir el modelo original para usar en visualizaciones
modelo_original <- lm(mpg ~ wt + hp, data = mtcars)

# =============================================================================
# 6. VISUALIZACIONES
# =============================================================================

# Crear gráficos
p1 <- ggplot(data.frame(medias = medias_bootstrap), aes(x = medias)) +
  geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7, color = "black") +
  geom_vline(xintercept = mean(mtcars$mpg), color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = quantile(medias_bootstrap, c(0.025, 0.975)), 
             color = "blue", linetype = "dotted", linewidth = 1) +
  labs(title = "Distribución Bootstrap de la Media MPG",
       x = "Media MPG", y = "Frecuencia") +
  theme_minimal()

# Gráfico Q-Q para normalidad
p2 <- ggplot(data.frame(medias = medias_bootstrap), aes(sample = medias)) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot: Distribución Bootstrap vs Normal",
       x = "Cuantiles Teóricos", y = "Cuantiles Muestrales") +
  theme_minimal()

# Bootstrap vs datos originales
datos_comparacion <- data.frame(
  valor = c(mtcars$mpg, medias_bootstrap),
  tipo = c(rep("Datos Originales", nrow(mtcars)), 
           rep("Medias Bootstrap", n_bootstrap))
)

p3 <- ggplot(datos_comparacion, aes(x = valor, fill = tipo)) +
  geom_density(alpha = 0.6) +
  labs(title = "Comparación: Datos Originales vs Distribución Bootstrap",
       x = "MPG", y = "Densidad") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "blue"))

# Convergencia del bootstrap
medias_acumuladas <- cumsum(medias_bootstrap) / seq_along(medias_bootstrap)
p4 <- ggplot(data.frame(iteracion = 1:n_bootstrap, media_acum = medias_acumuladas),
             aes(x = iteracion, y = media_acum)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = mean(mtcars$mpg), color = "red", linetype = "dashed") +
  labs(title = "Convergencia del Bootstrap",
       x = "Número de Iteraciones", y = "Media Acumulada") +
  theme_minimal()

# Mostrar todos los gráficos por separado para mejor visualización
print(p1)
print(p2) 
print(p3)
print(p4)

# También crear un gráfico combinado
tryCatch({
  grid.arrange(p1, p2, p3, p4, ncol = 2)
}, error = function(e) {
  cat("Error al crear gráfico combinado:", e$message, "\n")
  cat("Los gráficos se mostraron individualmente arriba.\n")
})

# =============================================================================
# 7. GRÁFICOS ADICIONALES MEJORADOS
# =============================================================================

# Gráfico de distribución de coeficientes de regresión
coefs_df <- data.frame(
  Intercepto = bootstrap_reg$t[,1],
  Peso_wt = bootstrap_reg$t[,2], 
  HP = bootstrap_reg$t[,3]
)

p5 <- ggplot(coefs_df, aes(x = Intercepto)) +
  geom_histogram(bins = 40, fill = "lightgreen", alpha = 0.7, color = "black") +
  geom_vline(xintercept = coef(modelo_original)[1], color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribución Bootstrap del Intercepto",
       x = "Valor del Intercepto", y = "Frecuencia") +
  theme_minimal()

p6 <- ggplot(coefs_df, aes(x = Peso_wt)) +
  geom_histogram(bins = 40, fill = "orange", alpha = 0.7, color = "black") +
  geom_vline(xintercept = coef(modelo_original)[2], color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribución Bootstrap del Coeficiente de Peso",
       x = "Coeficiente de Peso (wt)", y = "Frecuencia") +
  theme_minimal()

p7 <- ggplot(data.frame(correlaciones = bootstrap_cor$t), aes(x = correlaciones)) +
  geom_histogram(bins = 40, fill = "purple", alpha = 0.7, color = "black") +
  geom_vline(xintercept = cor(mtcars$mpg, mtcars$wt), color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribución Bootstrap de la Correlación MPG-WT",
       x = "Correlación", y = "Frecuencia") +
  theme_minimal()

# Scatterplot con líneas de regresión bootstrap
p8 <- ggplot(mtcars, aes(x = wt, y = mpg)) +
  geom_point(size = 3, alpha = 0.7, color = "darkblue") +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1.5) +
  labs(title = "Datos Originales: MPG vs Peso",
       x = "Peso (wt)", y = "Millas por Galón (mpg)") +
  theme_minimal()

print(p5)
print(p6)
print(p7)
print(p8)

# Definir el modelo original fuera del resumen para usarlo antes
# YA ESTÁ DEFINIDO ARRIBA

# =============================================================================
# 8. RESUMEN DE RESULTADOS
# =============================================================================

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("RESUMEN DE RESULTADOS BOOTSTRAP\n")
cat(paste(rep("=", 60), collapse=""), "\n")

cat("\nDATOS ORIGINALES:")
cat("\n- Tamaño de muestra:", nrow(mtcars))
cat("\n- Media MPG:", round(mean(mtcars$mpg), 3))
cat("\n- Desviación estándar MPG:", round(sd(mtcars$mpg), 3))
cat("\n- Correlación MPG-WT:", round(cor(mtcars$mpg, mtcars$wt), 3))

cat("\n\nRESULTADOS BOOTSTRAP (n =", n_bootstrap, "):")
cat("\n- Error estándar de la media:", round(sd(bootstrap_media$t), 3))
cat("\n- IC 95% para la media:", round(quantile(bootstrap_media$t, c(0.025, 0.975)), 3))
cat("\n- IC 95% para correlación MPG-WT:", round(ic_cor$percent[4:5], 3))

cat("\n\nCOEFICIENTE DE REGRESIÓN (mpg ~ wt + hp):")
cat("\n- Intercepto original:", round(coef(modelo_original)[1], 3))
cat("\n- IC 95% intercepto:", round(ic_intercepto$percent[4:5], 3))
cat("\n- Coef. wt original:", round(coef(modelo_original)[2], 3))
cat("\n- IC 95% wt:", round(ic_wt$percent[4:5], 3))

# =============================================================================
# 9. COMPARACIÓN DE MÉTODOS DE IC
# =============================================================================

cat("\n\nCOMPARACIÓN DE INTERVALOS DE CONFIANZA PARA LA MEDIA:")
cat("\n- Normal:", round(ic_media$normal[2:3], 3))
cat("\n- Básico:", round(ic_media$basic[4:5], 3))
cat("\n- Percentil:", round(ic_media$percent[4:5], 3))
cat("\n- BCa:", round(ic_media$bca[4:5], 3))

# Estadísticas adicionales
cat("\n\nESTADÍSTICAS ADICIONALES:")
cat("\n- Sesgo estimado (media):", round(bootstrap_media$t0 - mean(bootstrap_media$t), 4))
cat("\n- Error estándar teórico:", round(sd(mtcars$mpg)/sqrt(nrow(mtcars)), 3))
cat("\n- Error estándar bootstrap:", round(sd(bootstrap_media$t), 3))
cat("\n- Ratio EE_bootstrap/EE_teorico:", round(sd(bootstrap_media$t)/(sd(mtcars$mpg)/sqrt(nrow(mtcars))), 3))

cat("\n\n¡FIN!\n")