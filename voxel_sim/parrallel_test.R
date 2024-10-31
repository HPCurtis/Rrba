library("cmdstanr")
library("tidyverse")
library("posterior")

options(mc.cores = parallel::detectCores())  # Use multiple cores

# Specify file path. 
# TODO upload data to github as rba.csv
FILE_PATH <- "../data/df2.csv"
MOD_FILE_PATH <- "../cmdstanr_implementation/stanfiles/rba_parrallel.stan"

# Read in the fMRI Region BOLD data.
df <- read.csv(FILE_PATH)
print(nrow(df))

# Compile Stan model
mod <- cmdstan_model(MOD_FILE_PATH, cpp_options = list(stan_threads = TRUE), compile = TRUE)

# convert to subject and ROI to integer.
df$int_subj <- as.integer(factor(df$participant+1))
df$int_roi <- as.integer(factor(df$voxel_id +1))

X <- model.matrix(BOLD ~ condition, data= df)

# Names correspond to the data block in the Stan program
data_list <- list(N = nrow(df), y = df$BOLD, X = X, K = ncol(X), Kc = ncol(X) -1,
                  J = 2, N_subj = length(unique(df$participant)),
                  N_ROI=length(unique(df$voxel_id)), subj = df$int_subj, ROI=df$int_roi, 
                  # this value very fast but more divergences than n/8
                  grainsize=round(nrow(df)/16))

# Fit model.
fit <- mod$sample(
  data = data_list,
  seed = 123,
  iter_warmup = 1000,                  
  iter_sampling = 1000, 
  parallel_chains = 4,
  threads_per_chain = 2
)

# Extract posterior draws
posterior_df <- fit$draws(format = "df")

# Alternatively, to store them in a CSV file
write.csv(posterior_df, "posterior_parrallel_draws.csv", row.names = FALSE)

# Output summary to check for convergence.
# rhat is fine and ess_bulk & tail > 400.
print(fit$summary(variables = c("alpha", "b_Intercept", "tau_u")))
