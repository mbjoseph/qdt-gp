data {
  int<lower=1> n;
  matrix[n, n] D_sq;
  vector[n] z;
}
transformed data {
  matrix[n, n] jitter;
  jitter <- diag_matrix(rep_vector(0.00001, n));
}
parameters {
  real mu; 
  real<lower=0> eta_sq;
  real<lower=0> phi;
}
transformed parameters {
  cov_matrix[n] Sigma;
  Sigma <- eta_sq * exp(-phi * D_sq) + jitter;
}
model {
  mu ~ normal(0, 1);
  eta_sq ~ normal(0, 5);
  phi ~ normal(0, 5);
  z ~ multi_normal(rep_vector(mu, n), Sigma);
}