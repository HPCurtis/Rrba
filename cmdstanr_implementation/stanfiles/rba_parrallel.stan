functions {
  /* integer sequence of values
   * Args:
   *   start: starting integer
   *   end: ending integer
   * Returns:
   *   an integer sequence from start to end
   */
  array[] int sequence(int start, int end) {
    array[end - start + 1] int seq;
    for (n in 1:num_elements(seq)) {
      seq[n] = n + start - 1;
    }
    return seq;
  }
  // compute partial sums of the log-likelihood
  real partial_log_lik_lpmf(array[] int seq, int start, 
    int end, data vector y, data matrix X, data matrix Xc, vector beta, real alpha,
    real sigma, data array[] int subj, data array[] int ROI,
    vector u, matrix u2) {

    real ptarget = 0;
    int N = end - start + 1;
    // initialize linear predictor term
    vector[N] mu = alpha + u[subj[start:end]] + u2[ROI[start:end], 1] + X[,2][start:end] .* u2[ROI[start:end], 2];
    ptarget += normal_id_glm_lpdf(y[start:end] | Xc[start:end], mu, beta, sigma);
    return ptarget;
  }
}

data {
  int<lower=1> N;                            // Total number of observations
  int<lower=1> J;                            // Number of groups (subjects or ROIs)
  int<lower=1> N_subj;                       // Number of subjects
  int<lower=1> N_ROI;                        // Number of Regions of Interest (ROIs)
  int<lower=1> K;                            // Number of population-level predictors
  int<lower=1> Kc;                           // Number of centered population-level predictors
  int grainsize;
  vector[N] y;                               // Response variable
  matrix[N, K] X;                            // Design matrix for predictors
  array[N] int<lower=1, upper=N_subj> subj; // Subject IDs for each observation
  array[N] int<lower=1, upper=N_ROI> ROI;   // ROI IDs for each observation
}

transformed data {
  matrix[N, Kc] Xc;                          // Centered version of X without an intercept
  vector[Kc] means_X;                        // Column means of X before centering
  array[N] int seq = sequence(1, N);
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
  u2 = transpose((diag_pre_multiply(tau_u2, L_u) * z_u2));                // Centered random effects for ROIs

  target += reduce_sum(partial_log_lik_lpmf,seq, grainsize, y, X, Xc, 
                       beta, alpha, sigma, subj, ROI, u, u2);
}
