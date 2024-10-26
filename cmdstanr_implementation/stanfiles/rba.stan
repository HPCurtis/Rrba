data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower = 1> N_subj;
  int<lower = 1> K;
  int<lower = 1> Kc;

  vector[N] y;
  matrix[N, K] x;
  array[N] int<lower = 1, upper = N_subj> subj;
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
  vector<lower = 0>[J]  tau_u;
  real alpha;
  real beta;
  matrix[J, N_subj] z_u;
  cholesky_factor_corr[J] L_u;
  real<lower = 0> sigma;
}
transformed parameters {
  matrix[N_subj, J] u;
  u = (diag_pre_multiply(tau_u, L_u) * z_u)';
}
model {
  target += student_t_lpdf(alpha | 3, 0.1, 2.5);
  target += normal_lpdf(beta | 0, 10);
  target += student_t_lpdf(sigma | 3, 0, 2.5) 
  - student_t_lccdf(0 | 3, 0, 2.5);
  
  target += student_t_lpdf(tau_u | 3, 0, 2.5)
  -  student_t_lccdf(0 | 3, 0, 2.5);
  target += lkj_corr_cholesky_lpdf(L_u | 1);
  target += std_normal_lpdf(to_vector(z_u));

  // Generate model mu.
  vector[N] mu = alpha + u[subj, 1] + x .* (beta + u[subj, 2]);
  target += normal_lpdf(y | mu, sigma);
}
