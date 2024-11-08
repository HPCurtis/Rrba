# Load in analysis packages 
library("brms")
library("rstan")
library("tidyverse")
library("cmdstanr")

options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)

MOD_FILE_PATH = "brms_model_opt.stan"
FILE_PATH = "../data/data.csv"
df <- read.csv(FILE_PATH)

# Define your model formula
mod_formula <- bf(y ~ x + (1 | subject) + (x | ROI))

# Genrate the stan code and data for understanding.
sc <- stancode(mod_formula, data = df)
sc_parrallel <- stancode(mod_formula, data = df)
sd <- standata(mod_formula, data = df)
sd_p <- standata(mod_formula, data = df)

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

# Fit using AFNI defaults brms models (see, https://github.com/afni/afni/blob/master/src/R_scripts/RBA.R).
fit <- brm(
  formula = mod_formula,      
  data = df,
  init=0,
  control = list(adapt_delta = 0.99, max_treedepth = 15),                
  family = gaussian(),        
  backend = "cmdstanr",       
  chains = 4)

# Use brms multiple threading for speed ups.
fit_parrallel <- brm(
  formula = mod_formula,      
  data = df,
  init=0,
  control = list(adapt_delta = 0.99, max_treedepth = 15),                  
  family = gaussian(),        
  backend = "cmdstanr",       
  chains = 4,  
  cores = 4,               
  threads = threading(2))