# Simulating Gaussian processes
library(akima)
library(scales)

## Univariate inputs ---------------------------
n <- 1000
x <- sort(runif(n, -10, 10))

# squared exponential --------------------------
l <- 10  # range parameter
dmat <- as.matrix(dist(x))
C <- exp(-dmat^2 / 2 * l^2) + diag(rep(1E-6, n))

# z %*% L produces multivariate normal draws from MVN(0, Sigma), 
# where L %*% t(L) = Sigma. i.e., L is a cholesky decomposition of Sigma
# and z ~ normal(0, 1)
z <- rnorm(n)
y <- z %*% chol(C + d)
plot(x, y, type='l')

# or, built into a linear model
alpha <- -1
beta <- .5
y <- alpha + beta * x + z %*% chol(C + d)
plot(x, y, type='l')

# OU covariance function (more jaggety) --------------
C <- exp(-dmat / l)
z <- rnorm(n)
y <- z %*% chol(C + d)
plot(x, y, type='l')


## 2d gaussian process -------------------------------
# squared exponential covariance function:
C <- function(sigma, d, rho){
  stopifnot(sigma > 0)
  stopifnot(rho > 0)
  sigma ^ 2 * exp(-d^2 / (2 * rho ^ 2))
}

sigma <- 1
rho <- .7
sigma_e <- .001
n <- 2000
x1 <- runif(n, 0, 10)
x2 <- runif(n, 0, 10)
d_mat <- as.matrix(dist(cbind(x1, x2)))
Cmat <- C(sigma, d_mat, rho)
Emat <- diag(rep(sigma_e^2, n))
# simulate realizations
L_c <- t(chol(Cmat + Emat))
z <- rnorm(n)
y <- L_c %*% z
s <- interp(x1, x2, y, nx=200, ny=200)
image(s, col=rainbow(50))
points(x1, x2, col=alpha(1, .2), cex=.2)

# or an interactive 3d version
library(rgl)
persp3d(s$x, s$y, s$z, col = "lightblue")
