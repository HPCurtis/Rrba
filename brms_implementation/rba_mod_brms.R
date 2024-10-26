library("brms")
library("tidyverse")
options(mc.cores = parallel::detectCores())

file_path = "../data/data.csv"
df <- read.csv(file_path)
print(head(df))

# Define your model formula
mod_formula <- bf(y ~ x + (1 | subject) + (x | ROI))

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

sc <- stancode(fit)
sc_parrallel <- stancode(fit_parrallel)

# Define the file path for where to save brms outputted stan code.
stan_file_path <- "brms_model.stan"  
stan_file_path_parrallel <- "brms_model_parrallel.stan" 

# Write the Stan code to a .stan file
writeLines(sc, con = stan_file_path)
writeLines(sc_parrallel, con = stan_file_path_parrallel)