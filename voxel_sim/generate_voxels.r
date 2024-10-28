# Load necessary libraries
library(dplyr)
library(tidyr)

# Set the random seed for reproducibility
set.seed(42)

# Create voxel noise parameters
noise_list <- data.frame(
  mean = rnorm(1000, 0, 0.08),
  sd = abs(rnorm(1000, 0, 0.05))
)

# Function to generate random voxel values
generate_random_voxels <- function(mean, sd, noise_list, length = 1000) {
  voxels <- numeric(length)
  for (v in seq_len(length)) { # Adjusted for 1-based indexing
    mean_adj <- mean + noise_list$mean[v]
    sd_adj <- sd + noise_list$sd[v]
    voxels[v] <- rnorm(1, mean = mean_adj, sd = sd_adj)
  }
  return(voxels)
}

# Create participant and condition vectors for multi-level index
participants <- rep(1:30, each = 2) # Adjusted for 1-based indexing
conditions <- rep(1:2, times = 30) - 1  # Condition now matches Python's 0 and 1 values

# Initialize the data matrix
data_matrix <- matrix(0, nrow = 60, ncol = 1000)
df <- as.data.frame(data_matrix)

# Populate the data matrix with voxel values for each participant and condition
for (participant in 1:30) { # 1-based indexing for participants
  random_effect_mean <- rnorm(1, 0, 0.1)
  random_effect_sd <- abs(rnorm(1, 0, 0.05))
  
  for (condition in 1:2) { # 1-based indexing for conditions
    mean <- if (condition == 1) 0.5 + random_effect_mean else 0.3 + random_effect_mean
    sd <- 0.1 + random_effect_sd
    row_index <- (participant - 1) * 2 + condition # Calculate row index
    df[row_index, ] <- generate_random_voxels(mean, sd, noise_list)
  }
}

# Add participant and condition columns to match multi-level index
df <- df %>%
  mutate(participant = participants, condition = conditions) %>%
  pivot_longer(cols = -c(participant, condition), names_to = "voxel_id", values_to = "BOLD") %>%
  mutate(voxel_id = as.integer(sub("V", "", voxel_id)))

# Write the dataframe to a CSV file
write.csv(df, "df2.csv", row.names = FALSE)