# Install required packages if not installed
install.packages(c("reticulate", "ggplot2", "reshape2", "viridis"))

# Load necessary libraries
library(reticulate)
library(ggplot2)
library(reshape2)
library(viridis)

# Import NumPy and load the .npy file
np <- import("numpy")
data <- np$load("photon_decay_data.npy")  # Change path if needed

# Check dimensions of the array
dim(data)  # Should return (40, 60, 50) -> (rows, cols, time slices)

# Summary statistics of the entire dataset
summary(data)

# Extract a specific time frame (slice)
slice_1 <- data[,,1]  # Extract the first slice
dim(slice_1)  # Should be (40, 60), a 2D matrix

# Convert to a data frame for analysis
df <- as.data.frame(slice_1)

# Convert data to long format for ggplot visualization
df_long <- melt(slice_1)

# Heatmap of the first slice
ggplot(df_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Photon Decay Data - First Time Frame",
       x = "Pixel Column", y = "Pixel Row", fill = "Intensity")

# Iterate Over Slices for Summary Statistics
for (i in 1:dim(data)[3]) {
  cat("Slice:", i, "\n")
  print(summary(data[,,i]))
}





# Install required packages if not installed
install.packages(c("gganimate", "gifski", "png"))

# Load necessary libraries
library(ggplot2)
library(gganimate)
library(gifski)
library(png)
library(reshape2)
library(viridis)
library(reticulate)

# Import NumPy and load the .npy file
np <- import("numpy")
data <- np$load("photon_decay_data.npy")  # Adjust path if needed

# Convert dataset to long format for animation
df_long_time <- melt(data, varnames = c("Row", "Col", "Time"))

# Create the animation
p <- ggplot(df_long_time, aes(x = Row, y = Col, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  transition_time(Time) + 
  labs(title = "Photon Decay Over Time: Frame {frame}")

# Render and display the animation
animate(p, renderer = gifski_renderer(), fps = 10, width = 600, height = 400)

# Save animation as a GIF (optional)
anim_save("photon_decay_animation.gif", animation = p, renderer = gifski_renderer())

