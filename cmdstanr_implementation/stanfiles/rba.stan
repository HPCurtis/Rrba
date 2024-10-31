data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower = 1> N_subj;
  int<lower = 1> N_ROI;
  int<lower = 1> K;
  int<lower = 1> Kc;

  vector[N] y;
  matrix[N, K] X;
  array[N] int<lower = 1, upper = N_subj> subj;
  array[N] int<lower = 1, upper = N_ROI> ROI;
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  real<lower = 0> tau_u;
  vector<lower = 0>[J] tau_u2;
  real alpha;
  vector[Kc] beta;
  vector[N_subj] z_u;
  matrix[J, N_ROI] z_u2;
  cholesky_factor_corr[J] L_u;
  real<lower = 0> sigma;
}
model {
  target += student_t_lpdf(alpha | 3, 0.1, 2.5);
  target += normal_lpdf(beta | 0, 10);
  target += student_t_lpdf(sigma | 3, 0, 2.5) 
  - student_t_lccdf(0 | 3, 0, 2.5);
  
  target += student_t_lpdf(tau_u | 3, 0, 2.5)
  - student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(tau_u2 | 3, 0, 2.5)
  - student_t_lccdf(0 | 3, 0, 2.5);
  target += lkj_corr_cholesky_lpdf(L_u | 1);
  target += std_normal_lpdf(to_vector(z_u2));
  target += std_normal_lpdf(z_u);

  vector[N_subj] u;
  matrix[N_ROI, J] u2;
  u = z_u * tau_u;
  u2 = transpose((diag_pre_multiply(tau_u2, L_u) * z_u2));

  // Generate model mu.
  vector[N] mu = alpha + u[subj] + u2[ROI, 1] + X[,2] .* u2[ROI, 2];

  target += normal_id_glm_lpdf(y | Xc, mu, beta, sigma);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = alpha - dot_product(means_X, beta);
  
  // compute group-level correlations
  corr_matrix[J] Cor_1 = multiply_lower_tri_self_transpose(L_u);
  vector<lower=-1,upper=1>[1] cor_1;
   extract upper diagonal of correlation matrix
  for (k in 1:J) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
