# fitting a univariate gaussian process model
library(scales)
source('hdi.R')
par(bty='n')

# define some arbitrary function and simulate data with error
f <- function(x){
  # squiggly ones
  x * sin(x)
  #x - x^2 + -80 * exp(-x) + .8 * (x - mean(x))^3
  
  # less squiggly ones
  #.2 * x
  #-5*exp(-x)
}
n <- 20
x <- sort(runif(n, 0, 10))
error <- rnorm(n, 0, .2)
y <- f(x) + error 

plot(x, y, col='red')
lines(seq(min(x), max(x), length.out=200), 
      f(seq(min(x), max(x), length.out=200)))
y <- c(scale(y))

# simulate missingness in response vector
p_miss <- .4
n_miss <- round(p_miss * n)
obs_rm <- sort(sample(n, n_miss))
x_obs <- x[-obs_rm]
y_obs <- y[-obs_rm]
plot(x_obs, y_obs, 
     xlab='x', ylab='y', 
     ylim=c(-5, 5), xlim=c(-1, 11), pch=19)

# fit squared exponential GP model
library(rstan)
stan_d <- list(n = n - n_miss, 
               t = x_obs, 
               y = y_obs, 
               n_pred = n_miss, 
               t_pred = x[obs_rm])
watch <- c('eta_sq', 'phi', 'mu', 'sigma_sq', 'y_pred')
m_init <- stan('1d_sq_exp.stan', 
               data=stan_d, chains=1, iter=10)
m_fit <- stan(fit=m_init, chains=2, cores=2, 
              iter=1000, data=stan_d, pars=watch)
m_fit
post <- extract(m_fit)

# plot the HDIs for predictions estimated with Stan
y_pred_med <- apply(post$y_pred, 2, median)
n_iter <- length(post$lp__)
for (i in 1:n_miss){
  q <- quantile(post$y_pred[, i], probs=c(0.025, 0.975))
  segments(x0=x[obs_rm][i], x1=x[obs_rm][i], y0=q[1], y1=q[2], col='magenta')
}

# Interpolation & the posterior predictive density ----------------------------
# e_y calculates the conditional expectation of new data
e_y <- function(mu_1, omega_12, omega_22, y_2, mu_2){
  mu_1 + omega_12 %*% solve(omega_22) %*% (y_2 - mu_2)
}
# var_y calculates the conditional variance of new data
var_y <- function(omega_11, omega_12, omega_22, omega_21){
  omega_11 - omega_12 %*% solve(omega_22) %*% omega_21
}

new_x <- seq(-1, 11, .05)
n_new <- length(new_x)
all_x <- c(new_x, x_obs)
nt <- length(all_x)
D <- as.matrix(dist(all_x)^2)
y_2 <- y_obs
e_new <- array(dim=c(n_new, n_iter))
var_new <- array(dim=c(n_new, n_iter))

# cycle through iterations and calculate posterior E(x_new) & Var(x_new)
pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
for (i in 1:n_iter){
  eta_sq <- post$eta_sq[i]
  phi <- post$phi[i]
  Omega <- eta_sq * exp(-phi * D) + diag(rep(post$sigma_sq[i], nt))
  omega_11 <- Omega[1:n_new, 1:n_new]
  omega_22 <- Omega[(n_new+1):nt, (n_new+1):nt]
  omega_12 <- Omega[1:n_new, (n_new+1):nt]
  omega_21 <- t(omega_12)
  mu_1 <- post$mu[i]
  mu_2 <- post$mu[i]
  e_new[, i] <- e_y(mu_1, omega_12, omega_22, y_2, mu_2)
  var_new[, i] <- diag(var_y(omega_11, omega_12, omega_22, omega_21))
  setTxtProgressBar(pb, i)
}

# draw posterior expected values
mu_e <- apply(e_new, 1, mean)
hdi_e <- apply(e_new, 1, hdi)
lines(new_x, mu_e, col='blue')
lines(new_x, hdi_e[1, ], col='blue', lty=2)
lines(new_x, hdi_e[2, ], col='blue', lty=2)
par(xpd=TRUE)
legend(2, 4, lty=c(1, 1), col=c('blue', 'magenta'), 
       legend=c(expression(paste('P(', mu[y], '|', y[obs], ')')),
                expression(paste('P(', tilde(y), '|', y[obs], ')'))), 
       bty='n', inset=c(0,-.3), horiz=TRUE)
par(xpd=FALSE)

# posterior predictive distribution
ntotal <- length(c(e_new))
# marginals for y_new are univariate normal
sd_new <- sqrt(var_new)
y_new <- array(rnorm(ntotal, e_new, sd_new), dim=dim(e_new))
y_hdi <- apply(y_new, 1, hdi)
lines(new_x, y_hdi[1, ], col='magenta')
lines(new_x, y_hdi[2, ], col='magenta')
#points(x[obs_rm], y[obs_rm], pch=19, col='green') # uncomment: true missing vals
