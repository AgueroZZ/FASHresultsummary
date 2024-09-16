# write a function to simulate process given sd_fun and boundary_sd
simulate_data <- function(x, a, k = 30, sd_fun, boundary_sd, intercept_sd = 0.5, n){
  m <- 1
  P <- BayesGP:::Compute_Q_sB(a = a, k = k, region = c(0,2), accuracy = 1000)
  B <- BayesGP:::Compute_B_sB_helper(refined_x = x, a = a, k = k, m = m, region = c(0,2))

  ## simulate from precision matrix
  random_weights <- LaplacesDemon::rmvnp(1, mu = rep(0, ncol(P)), Omega = as.matrix(((1/sd_fun)^2) * P))
  f_standard <- B %*% as.vector(random_weights)

  ## add boundary condition
  boundary_weights <- rnorm(m, 0, boundary_sd)
  f <- f_standard + sin(a*x) * boundary_weights

  # add intercept
  intercept <- rnorm(1, 0, intercept_sd)

  y <- rpois(lambda = exp(intercept + f), n = n)
  return(list(x = x, y = y, f = f))
}

# second function to simulate quasi-periodic data
simulate_data_v2 <- function(x, a, sd_quasi = 1, sd_exact = 1, intercept_sd = 0.5, n){
  f_standard <- sin(a*x) * rnorm(1, 0, sd_exact) + cos(a*x) * rnorm(1, 0, sd_exact)
  f_quasi <- sin(a*x) * x * rnorm(1, 0, sd_quasi) + cos(a*x) * x * rnorm(1, 0, sd_quasi)
  f <- f_standard + f_quasi
  # add intercept
  intercept <- rnorm(1, 0, intercept_sd)

  y <- rpois(lambda = exp(intercept + f), n = n)
  return(list(x = x, y = y, f = f))
}
