functions {
  // Function to compute log likelihood for a subset of data
  real log_likelihood(int start, int end,
                      vector y, matrix X, 
                      vector alpha, vector beta,
                      vector u, matrix u2,
                      real sigma) {
    real ll = 0;
    
    for (n in start:end) {
      vector[Kc] mu; 
      mu = alpha + u[subj[n]] + u2[ROI[n], 1] + X[n, 2] * u2[ROI[n], 2];
      ll += normal_id_glm_lpdf(y[n] | mu, beta, sigma);
    }
    return ll;
  }
}

data {
  int<lower=1> N;                            // Total number of observations
  int<lower=1> J;                            // Number of groups (subjects or ROIs)
  int<lower=1> N_subj;                       // Number of subjects
  int<lower=1> N_ROI;                        // Number of Regions of Interest (ROIs)
  int<lower=1> K;                            // Number of population-level predictors
  int<lower=1> Kc;                           // Number of centered population-level predictors

  vector[N] y;                               // Response variable
  matrix[N, K] X;                            // Design matrix for predictors
  array[N] int<lower=1, upper=N_subj> subj; // Subject IDs for each observation
  array[N] int<lower=1, upper=N_ROI> ROI;   // ROI IDs for each observation
}

transformed data {
  matrix[N, Kc] Xc;                          // Centered version of X without an intercept
  vector[Kc] means_X;                        // Column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}

parameters {
  real<lower = 0> tau_u;                     // Standard deviation for subject-level effects
  vector<lower = 0>[J] tau_u2;               // Standard deviations for ROI effects
  real alpha;                                 // Population-level intercept
  vector[Kc] beta;                            // Coefficients for the centered predictors
  vector[N_subj] z_u;                         // Standardized group-level effects for subjects
  matrix[J, N_ROI] z_u2;                      // Standardized group-level effects for ROIs
  cholesky_factor_corr[J] L_u;                // Cholesky factor of the correlation matrix for u2
  real<lower = 0> sigma;                      // Residual standard deviation
}

model {
  target += student_t_lpdf(alpha | 3, 0.1, 2.5);                // Prior for alpha
  target += normal_lpdf(beta | 0, 10);                           // Prior for beta
  target += student_t_lpdf(sigma | 3, 0, 2.5) 
            - student_t_lccdf(0 | 3, 0, 2.5);                    // Prior for sigma
  target += student_t_lpdf(tau_u | 3, 0, 2.5) 
            - student_t_lccdf(0 | 3, 0, 2.5);                    // Prior for tau_u
  target += student_t_lpdf(tau_u2 | 3, 0, 2.5) 
            - student_t_lccdf(0 | 3, 0, 2.5);                    // Prior for tau_u2
  target += lkj_corr_cholesky_lpdf(L_u | 1);                     // Prior for correlation matrix
  target += std_normal_lpdf(to_vector(z_u2));                    // Prior for group-level effects of ROIs
  target += std_normal_lpdf(z_u);                                 // Prior for group-level effects of subjects

  // Compute random effects
  vector[N_subj] u;
  matrix[N_ROI, J] u2;
  u = z_u * tau_u;                                                // Centered random effects for subjects
  u2 = (diag_pre_multiply(tau_u2, L_u) * z_u2)';                // Centered random effects for ROIs

  // Compute the log likelihood using reduce_sum
  target += reduce_sum(log_likelihood, 
                        to_array_1d(y), 
                        to_array_1d(X), 
                        alpha, beta, u, u2, sigma,
                        1, // This is the chunk size, adjust based on memory limits
                        1, N); // Start and end indices for the likelihood
}