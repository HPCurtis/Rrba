library("cmdstanr")
library("tidyverse")

options(mc.cores = parallel::detectCores())  # Use multiple cores

# Specify file path. 
# TODO upload data to github as rba.csv
FILE_PATH <- "https://raw.githubusercontent.com/HPCurtis/Datasets/refs/heads/main/rba.csv"
MOD_FILE_PATH <- "stanfiles/rba.stan"

# Read in the fMRI Region BOLD data.
df <- read.csv(FILE_PATH)

# Compile stan model
mod <- cmdstan_model(MOD_FILE_PATH, cpp_options = list(stan_threads = TRUE), compile = TRUE)

# convert to subject and ROI to integer.
df$int_subj <- as.integer(factor(df$subject))
df$int_roi <- as.integer(factor(df$ROI))

X <- model.matrix(y ~ x, data= df)

# Names correspond to the data block in the Stan program
data_list <- list(N = nrow(df), y = df$y, X = X, K = ncol(X), Kc = ncol(X) -1,
                  J = 2, N_subj = length(unique(df$subject)),
                  N_ROI=length(unique(df$ROI)), subj = df$int_subj, ROI=df$int_roi)

# Fit model.
fit <- mod$sample(
  data = data_list,
  iter_warmup = 1000,                  
  iter_sampling = 1000,
  chains=4, 
  parallel_chains = 4,
  threads_per_chain = 2,
  adapt_delta=.99, 
  max_treedepth=15
)

fit$save_object(file = "RBA_posterior_draws.RDS")

# Extract posterior draws
posterior_df <- fit$draws(format = "df")

# Output summary to check for convergence.
# rhat is fine and ess_bulk & tail > 400.
print(fit$summary(variables = c("alpha", "beta", "sigma", "tau_u", "tau_u2")))
