

sim <- function(n.ind = 100, n.measure = 10, 
                x.d, coef.d = 1, tau.d = 1.5, 
                x.z, coef.z = 0, tau.z = 0, p.zero = 0,
                x.w, coef.w = 0, theta = 2)
{
  n <- n.ind * n.measure                                  # total number of observations
  ind <- sort(rep(1:n.ind, n.measure))                    # subject ID
  if (missing(x.d)) x.d <- c(rep(0, n/2), rep(1, n/2))    # a binary variable: disease vs. healthy
  x.d <- as.matrix(x.d)
  if (missing(x.z)) x.z <- x.d        
  x.z <- as.matrix(x.z)
  if (missing(x.w)) x.w <- x.d 
  x.w <- as.matrix(x.w)
  
  mu <- runif(n, 0.1, 3.5)
  tn <- exp(mu)/exp(-7)    

  b <- rep(NA, n.ind)                  # random effect: b ~ N(0, tau.z^2)
  for (j in 1:n.ind) b[j] <- rnorm(1, 0, tau.z)
  eta <- b[ind] + x.z %*% coef.z
  y.normal <- rnorm(n, -eta, 1.6)
  quantiles <- quantile(y.normal, p.zero)
  y.z <- as.numeric( factor(cut(y.normal, breaks = c(-Inf, quantiles, Inf))) ) - 1 

  mu.theta <- log(theta)
  theta <- exp(mu.theta + x.w %*% coef.w)
  
  b <- rep(NA, n.ind)                  # random effect: b ~ N(0, tau.c^2)
  for (j in 1:n.ind) b[j] <- rnorm(1, 0, tau.d)
  eta <- mu + b[ind] + x.d %*% coef.d
  y.nb <- rnbinom(n, mu = exp(eta), size = theta)
  
  y <- ifelse(y.z == 0, 0, y.nb)
  
  list(ind.ID=ind, x.d=x.d, x.z=x.z, x.w=x.w, T=tn, y=y, theta=theta, y.z=y.z)
}
