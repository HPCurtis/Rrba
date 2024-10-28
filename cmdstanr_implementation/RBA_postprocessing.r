library(tidyverse)
source("postprocessing_utilities.r")

FILE_PATH <- "posterior_draws.csv"
DATA_PATH <- "../data/data.csv"

posterior_draws <- read_csv(FILE_PATH)
df <- read.csv(DATA_PATH)


cor_draws <- L_u_extract(posterior_draws)

tau_extract <- tau_extract(posterior_draws)

rho <- LLt_conv(posterior_draws)
beta <- standardise_inverse(posterior_draws, df$x)
