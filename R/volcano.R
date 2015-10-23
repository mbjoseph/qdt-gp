# interpolating the volcano data with a Gaussian process model
library(fields)
library(reshape2)
library(rstan)

# subsample/prep the volcano data -----------------------
p_miss <- .99
vd <- melt(volcano)
n_missing <- round(p_miss * nrow(vd))
which_keep <- sample(nrow(vd), nrow(vd) - n_missing)
obs <- vd[which_keep, ]
obs$value <- c(scale(obs$value))

# calculate squared distances between points
D <- as.matrix(dist(obs[c('Var1', 'Var2')]))
D <- D / sd(D)
D_sq <- D^2
stan_d <- list(n=nrow(obs), D_sq=D_sq, z=obs$value)

# fit model in stan
sm_init <- stan('sq_exp.stan', data=stan_d, chains=1, iter=1)
watch <- c('mu', 'eta_sq', 'phi')
m_fit <- stan(fit=sm_init, data=stan_d, chains=2, cores=2, iter=400, pars=watch)
m_fit
post <- rstan::extract(m_fit)

# Interpolation & the posterior predictive density ----------------------------
# e_y calculates the conditional expectation of new data
e_y <- function(mu_1, omega_12, omega_22, y_2, mu_2){
  mu_1 + omega_12 %*% solve(omega_22) %*% (y_2 - mu_2)
}

new_x1 <- seq(1, nrow(volcano), 1)
new_x2 <- seq(1, ncol(volcano), 1)
new_x <- expand.grid(Var1 = new_x1, Var2 = new_x2)
n_new <- nrow(new_x)
all_x <- rbind(new_x, obs[c('Var1', 'Var2')])
nt <- nrow(all_x)
D_all <- as.matrix(dist(all_x))
D_all <- D_all / sd(D_all)
D_all <- D_all^2
y_2 <- obs$value
n_iter <- length(post$lp__)
e_new <- array(dim=c(n_new, n_iter))

# cycle through iterations and calculate posterior E(x_new) & Var(x_new)
pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
for (i in 1:n_iter){
  eta_sq <- post$eta_sq[i]
  phi <- post$phi[i]
  Omega <- eta_sq * exp(-phi * D_all)
  omega_11 <- Omega[1:n_new, 1:n_new]
  omega_22 <- Omega[(n_new+1):nt, (n_new+1):nt]
  omega_12 <- Omega[1:n_new, (n_new+1):nt]
  omega_21 <- t(omega_12)
  mu_1 <- post$mu[i]
  mu_2 <- post$mu[i]
  e_new[, i] <- e_y(mu_1, omega_12, omega_22, y_2, mu_2)
  setTxtProgressBar(pb, i)
}

# draw posterior expected values
mu_e <- apply(e_new, 1, mean)
hdi_e <- apply(e_new, 1, hdi, .99)
new_x$e_z <- mu_e
new_x$lo <- hdi_e[1, ]
new_x$hi <- hdi_e[2, ]

cd <- acast(new_x[c('Var1', 'Var2', 'e_z')], Var1 ~ Var2)
clo <- acast(new_x[c('Var1', 'Var2', 'lo')], Var1 ~ Var2)
chi <- acast(new_x[c('Var1', 'Var2', 'hi')], Var1 ~ Var2)

# evaluate predictive error
par(mfrow=c(2, 3))
par(mar=c(0, 1, 1, 1))
persp(volcano, theta=40, phi=45)
title('Volcano data', line=-1)
persp(cd, theta=40, phi=45)
title('Posterior mean based on 1% of data', line=-1)
persp((scale(volcano) - cd)^2, theta=40, phi=45)
title('Squared error', line=-1)
par(mar=c(5, 4, 1, 2) + 0.1)
image(volcano, col=terrain.colors(50))
points(obs$Var1 / nrow(volcano), obs$Var2 / ncol(volcano), pch=3)
image(cd, col=terrain.colors(50))
image.plot((scale(volcano) - cd)^2, col=topo.colors(50))
points(obs$Var1 / nrow(volcano), obs$Var2 / ncol(volcano), pch=3, col='white')
