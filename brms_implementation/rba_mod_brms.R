# Import analysis packages 
library("brms")
library("tidyverse")
options(mc.cores = parallel::detectCores())

file_path = "../data/data.csv"
df <- read.csv(file_path)

# Define your model formula
mod_formula <- bf(y ~ x + (1 | subject) + (x | ROI))

# Fit using brms defaults.
fit <- brm(
  formula = mod_formula,      
  data = df,                  
  family = gaussian(),        
  backend = "cmdstanr",       
  chains = 4,                 
  cores = 4)

# Use brms multiple threading for speed ups.
fit_parrallel <- brm(
  formula = mod_formula,      
  data = df,                  
  family = gaussian(),        
  backend = "cmdstanr",       
  chains = 4,                 
  cores = 4,
  threads = threading(2)      
)

# Gnerate the stan code and data for understanding.
sc <- stancode(fit)
sc_parrallel <- stancode(fit_parrallel)
sd <- standata(fit)
sd_p <- standata(fit_parrallel)

# Define the file path for where to save brms outputted stan code.
stan_file_path <- "brms_model.stan"  
stan_file_path_parrallel <- "brms_model_parrallel.stan" 
standata_file_path <- "standata.rds"
standata_file_path_parralel <- "standataparrallel.rds"

# Write the Stan code to a .stan file
writeLines(sc, con = stan_file_path)
writeLines(sc_parrallel, con = stan_file_path_parrallel)

# Create file for standata  
saveRDS(sd, standata_file_path)
saveRDS(sd_p, standata_file_path_parralel)