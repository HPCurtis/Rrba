library("cmdstanr")
library("tidyverse")

options(mc.cores = parallel::detectCores())  # Use multiple cores

# Specify file path. 
# TODO upload data to github as rba.csv
FILE_PATH <- "../data/data.csv"
MOD_FILE_PATH <- "stanfiles/rba.stan"

df <- read.csv(FILE_PATH)

# Compile stan model
mod <- cmdstan_model(MOD_FILE_PATH)

# convert to subject to integer.
df$int_subj <- as.integer(factor(df$subject))
df$int_roi <- as.integer(factor(df$ROI))

X <- model.matrix(y ~ x, data= df)


# Names correspond to the data block in the Stan program
data_list <- list(N = nrow(df), y = df$y, X = X, K = ncol(X), Kc = ncol(X) -1,
                  J = 2, N_subj = length(unique(df$subject)),
                  N_ROI=length(unique(df$ROI)), subj = df$int_subj, ROI=df$int_roi)

fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)