# Instalar y cargar librerías necesarias
install.packages(c("ggplot2", "reshape2", "viridis", "reticulate"))

library(ggplot2)
library(reshape2)
library(viridis)
library(reticulate)

# Importar NumPy y cargar el archivo .npy
np <- import("numpy")
data <- np$load("photon_decay_data.npy")  # Ajusta la ruta si es necesario

# Obtener dimensiones del conjunto de datos (Tsteps, Rows, Cols)
dim_data <- dim(data)
Tsteps <- dim_data[1]  # Número de frames (tiempo)
Rows <- dim_data[2]    # Número de filas (altura de la imagen)
Cols <- dim_data[3]    # Número de columnas (ancho de la imagen)

# Inicializar una matriz para almacenar las tasas de decaimiento estimadas
q_estimates <- matrix(NA, nrow = Rows, ncol = Cols)

# Bucle sobre cada píxel para estimar la tasa de decaimiento q
for (i in 1:Rows) {
  for (j in 1:Cols) {
    # Extraer la serie temporal del píxel (i, j)
    N_obs <- data[, i, j]  # Ahora el tiempo es la primera dimensión
    
    # Asegurar monotonía decreciente con el mínimo acumulado
    N_obs_adj <- cummin(N_obs)
    
    # Depuración: imprimir series temporales para verificar
    cat(sprintf("Pixel (%d, %d) - Original N_obs:\n", i, j), N_obs, "\n")
    cat(sprintf("Pixel (%d, %d) - Adjusted N_obs (cummin):\n", i, j), N_obs_adj, "\n")
    
    # Si la serie ajustada es completamente cero, asignar tasa de decaimiento cero
    if (all(N_obs_adj == 0)) {
      q_estimates[i, j] <- 0
      cat(sprintf("Pixel (%d, %d) - Solo ceros detectados. q = 0.\n\n", i, j))
    } else {
      # Calcular la tasa de decaimiento usando la fórmula MLE de un proceso exponencial
      num <- N_obs_adj[1] - N_obs_adj[Tsteps]  # Diferencia inicial-final
      den <- sum(N_obs_adj[1:(Tsteps - 1)])    # Suma de valores a riesgo
      q_hat <- ifelse(den < 1e-12, 0, num / den)  # Evitar divisiones por cero
      q_estimates[i, j] <- q_hat
      cat(sprintf("Pixel (%d, %d) - Estimated q_hat: %f\n\n", i, j, q_hat))
    }
  }
}

# Convertir la matriz q_estimates a formato largo para graficar con ggplot
df_q <- melt(q_estimates, varnames = c("Row", "Col"), value.name = "q_est")

# Crear un heatmap de las tasas de decaimiento estimadas
p_q <- ggplot(df_q, aes(x = Col, y = Row, fill = q_est)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Decay Rate q") +
  theme_minimal() +
  labs(title = "Estimated Decay Rates (Single Dye Model)")

# Mostrar la gráfica
print(p_q)
