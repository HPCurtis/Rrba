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
  matrix[N_ROI, J] z_u2;
  real<lower = 0> sigma;
}
transformed parameters{
  // Compute random effects
  vector[N_subj] u;
  matrix[N_ROI, J] u2;
  u = z_u * tau_u;                              
  u2[, 1] = z_u2[, 1] * tau_u2[1];
  u2[, 2] = z_u2[, 2] * tau_u2[2]; 
}
model {
  alpha ~ student_t(3, 0.1, 2.5);
  beta ~ normal(0, 1);
  sigma ~ student_t(3, 0, 2.5);
  
  tau_u ~ student_t(3, 0, 2.5);
  tau_u2 ~ student_t(3, 0, 2.5);
  
  to_vector(z_u2) ~ std_normal();
  z_u ~ std_normal();

  // Generate model mu.
  vector[N] mu = alpha + u[subj] + u2[ROI, 1] + X[,2] .* u2[ROI, 2];

  target += normal_id_glm_lpdf(y | Xc, mu, beta, sigma);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = alpha - dot_product(means_X, beta);
}
