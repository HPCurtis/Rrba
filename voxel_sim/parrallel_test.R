library("cmdstanr")
library("tidyverse")

options(mc.cores = parallel::detectCores())  # Use multiple cores

# Specify file path. 
# TODO upload data to github as rba.csv
FILE_PATH <- "df2.csv"
MOD_FILE_PATH <- "../cmdstanr_implmentation/stanfiles/rba_parrallel.stan"

# Read in the fMRI Region BOLD data.
df <- read.csv(FILE_PATH)
print(nrow(df))

# Compile stan model
mod <- cmdstan_model(MOD_FILE_PATH, cpp_options = list(stan_threads = TRUE))

# convert to subject and ROI to integer.
df$int_subj <- as.integer(factor(df$subject))
df$int_roi <- as.integer(factor(df$ROI))

X <- model.matrix(y ~ x, data= df)

# Names correspond to the data block in the Stan program
data_list <- list(N = nrow(df), y = df$y, X = X, K = ncol(X), Kc = ncol(X) -1,
                  J = 2, N_subj = length(unique(df$subject)),
                  N_ROI=length(unique(df$ROI)), subj = df$int_subj, ROI=df$int_roi, 
                  # brms default is N/4
                  grainsize=round(nrow(df)/4))

# Fit model.
fit <- mod$sample(
  data = data_list,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  threads_per_chain = 2
)