library(tidyverse)

L_u_extract <- function(posterior_df) {

    cor <- posteior_draws %>% 
            select("L_u[1,1]", "L_u[2,1]", "L_u[1,2]", "L_u[2,2]")
    return(cor)
}

tau_extract <- function(posterior_draws) {

    tau <- posterior_draws %>% 
            select("tau_u","tau_u2[1]","tau_u2[2]")
    return(tau)
}

LLt_conv <- function(posterior_draws) {
    cor_draws <- L_u_extract(posterior_draws)
    Cor_1 <- list()
    rho <- list()

for (i in 1:nrow(cor_draws)) {
    # Create the L_u matrix using the original column names
    L_u <- matrix(as.numeric(cor_draws[i, ]),
                  nrow = 2, ncol = 2, byrow = FALSE)
    # Compute LLt
    LLt <- L_u %*% t(L_u)
    Cor_1[[i]] <- LLt
    rho[[i]] <- Cor_1[[i]][upper.tri(Cor_1[[i]])]
}
    return(rho)
}

standardise_inverse <- function(posteriorr_draws, x) {
    beta <- posterior_draws %>% 
            select("beta[1]")
    beta <- (beta * sd(x)) + mean(x)
    return(beta)
}

