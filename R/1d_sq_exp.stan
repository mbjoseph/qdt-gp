data {
  int<lower=1> n;
  vector[n] t;
  vector[n] y;
  int<lower=1> n_pred;
  vector[n_pred] t_pred;
}
transformed data {
  int<lower=1> n_tot;
  n_tot <- n + n_pred;
}
parameters {
  real mu; 
  real<lower=0> eta_sq;
  real<lower=0> phi;
  real<lower=0> sigma_sq;
  vector[n_pred] y_pred;
}
transformed parameters {
  cov_matrix[n + n_pred] Sigma;
  vector[n_tot] t_all;
  
  for (i in 1:n) t_all[i] <- t[i];
  for (i in 1:n_pred) t_all[n + i] <- t_pred[i];
  
  // off-diagonal elements
  for (i in 1:(n_tot - 1)) {
    for (j in (i + 1):n_tot) {
      Sigma[i, j] <- eta_sq * 
                      exp(-phi * pow(t_all[i] - t_all[j], 2));
      Sigma[j, i] <- Sigma[i, j];
    }
  }
  
  // diagonal elements
  for (k in 1:n_tot)
    Sigma[k, k] <- eta_sq + sigma_sq; // + jitter for pos def.ness
}
model {
  vector[n_tot] y_all;

  mu ~ normal(0, 10);
  eta_sq ~ normal(0, 5);
  phi ~ normal(0, 5);
  sigma_sq ~ normal(0, 2);
 
  for (i in 1:n) y_all[i] <- y[i];
  for (i in 1:n_pred) y_all[i + n] <- y_pred[i];
  y_all ~ multi_normal(rep_vector(mu, n_tot), Sigma);
}