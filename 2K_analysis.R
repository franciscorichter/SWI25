# Install and load necessary libraries
install.packages(c("ggplot2", "reshape2", "gganimate", "viridis", "gifski", "png", "reticulate"))
library(ggplot2)
library(reshape2)
library(gganimate)
library(viridis)
library(gifski)
library(png)
library(reticulate)
library(Rcpp)
sourceCpp("mcem_multi_dye.cpp")
# Import NumPy and load the .npy file
np <- import("numpy")
data <- np$load("photon_decay_data.npy")  # Adjust path if needed

# Get dimensions of the dataset
dim_data <- dim(data)  # Should be (Rows, Cols, Time)
Rows <- dim_data[1]
Cols <- dim_data[2]
Tsteps <- dim_data[3]

# Parameters for MCEM
K <- 2  # Number of decay components (adjust as needed)
M <- 50  # Number of Monte Carlo samples
max_iter <- 100  # Maximum iterations
tol <- 1e-6  # Convergence tolerance
q_init <- rep(0.05, K)  # Initial values for q

# Initialize matrix to store estimated decay rates
q_estimates <- array(NA, dim = c(Rows, Cols, K))

# Loop over each pixel and apply MCEM
for (i in 1:Rows) {
  for (j in 1:Cols) {
    N_obs <- data[i, j, ]  # Extract time series for this pixel
    N_obs <- sort(N_obs,decreasing = T)
    # Print original time series for debugging
    cat(sprintf("Pixel (%d, %d) - Original N_obs:\n", i, j), N_obs, "\n")
    
    # Identify nonzero indices
    nonzero_idx <- which(N_obs > 0)
    
    if (length(nonzero_idx) > 0) {
      # Determine last nonzero frame and trim time series
      Tsteps_trim <- max(nonzero_idx)
      N_obs_trim <- N_obs[1:Tsteps_trim]
    } else {
      # If all zeros, retain original time steps
      Tsteps_trim <- Tsteps
      N_obs_trim <- N_obs
    }
    
    # Print trimmed time series for debugging
    cat(sprintf("Pixel (%d, %d) - Trimmed N_obs:\n", i, j), N_obs_trim, "\n")
    
    # Check if the trimmed dataset is entirely zero
    if (all(N_obs_trim == 0)) {
      cat(sprintf("Pixel (%d, %d) - Only zeros detected, skipping MCEM.\n", i, j))
      q_estimates[i, j, ] <- 0  # Assign zero decay rate if no signal
    } else {
      # Apply MCEM function
      est <- mcem_multi_dye_cpp(N_obs_trim, K, Tsteps_trim, q_init, M, max_iter, tol)
      
      # Print estimated q values
      cat(sprintf("Pixel (%d, %d) - Estimated q_hat:\n", i, j), est$q_hat, "\n")
      
      # Store estimated q values
      if(sum(sum(is.na(est$q_hat)))==0) q_estimates[i, j, ] <- sort(est$q_hat)
    }
  }
}


# Convert estimates to long format for visualization
df_q_long <- melt(q_estimates, varnames = c("Row", "Col", "Dye"), value.name = "q_est")

# Visualize estimated decay rates for the first dye
p_q <- ggplot(subset(df_q_long, Dye == 1), aes(x = Row, y = Col, fill = q_est)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Estimated Decay Rates for Dye 1", fill = "Decay Rate")

# Show heatmap of estimated decay rates
print(p_q)



p_q2 <- ggplot(subset(df_q_long, Dye == 2), aes(x = Row, y = Col, fill = q_est)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Estimated Decay Rates for Dye 2", fill = "Decay Rate")

# Show heatmap of estimated decay rates
print(p_q2)




# Compute difference between Dye 1 and Dye 2 decay rates
df_diff <- merge(
  subset(df_q_long, Dye == 1, select = c("Row", "Col", "q_est")),
  subset(df_q_long, Dye == 2, select = c("Row", "Col", "q_est")),
  by = c("Row", "Col"),
  suffixes = c("_Dye1", "_Dye2")
)

df_diff$q_diff <- df_diff$q_est_Dye1 - df_diff$q_est_Dye2  # Compute difference

# Create a heatmap for the difference
p_diff <- ggplot(df_diff, aes(x = Row, y = Col, fill = q_diff)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_minimal() +
  labs(title = "Difference in Decay Rates (Dye 1 - Dye 2)", fill = "Δ Decay Rate")

# Show the difference heatmap
print(p_diff)










