# write a function to simulate process given sd_fun and boundary_sd
simulate_data_v1 <- function(x, a, k = 30, sd_fun, boundary_sd, intercept_sd = 0.5){
  n <- length(x)
  m <- 1
  P <- BayesGP:::Compute_Q_sB(a = a, k = k, region = range(x), accuracy = 1000)
  B <- BayesGP:::Compute_B_sB_helper(refined_x = x, a = a, k = k, m = m, region = range(x))
  ## add boundary condition
  boundary_weights <- rnorm(m, 0, boundary_sd)
  # add intercept
  intercept <- rnorm(1, 0, intercept_sd)

  if(sd_fun != 0){
    sd_fun <- sd_fun/BayesGP:::compute_d_step_sgpsd(d = 1, a = a)
    ## simulate from precision matrix
    random_weights <- LaplacesDemon::rmvnp(1, mu = rep(0, ncol(P)), Omega = as.matrix(((1/sd_fun)^2) * P))
    f_standard <- B %*% as.vector(random_weights)

    f <- intercept + f_standard + sin(a*x) * boundary_weights
  } else{
    f <- intercept + sin(a*x) * boundary_weights
  }

  f <- f - mean(f) + 3

  y <- rpois(lambda = exp(f), n = n)
  return(list(x = x, y = y, f = f))
}

# second function to simulate quasi-periodic data
simulate_data_v2 <- function(x, a, sd_quasi = 1, sd_exact = 1, intercept_sd = 0.5){
  n <- length(x)
  f_standard <- sin(a*x) * rnorm(1, 0, sd_exact) + cos(a*x) * rnorm(1, 0, sd_exact)
  f_quasi <- sin(a*x) * x * rnorm(1, 0, sd_quasi) + cos(a*x) * x * rnorm(1, 0, sd_quasi)
  # add intercept
  intercept <- rnorm(1, 0, intercept_sd)

  f <- f_standard + f_quasi + intercept

  # for comparison purpose, center the mean at some value
  f <- f - mean(f) + 3

  y <- rpois(lambda = exp(f), n = n)
  return(list(x = x, y = y, f = f))
}

# simulate by choice

# simulate by choice
simulate_data <- function(type, x, a, k = 30, sd_fun = NULL, boundary_sd = NULL, intercept_sd = 0.5, sd_quasi = 1, sd_exact = 1, sd_np = 1, slope_sd = 0.5){

  if(type == "v1"){
    # For v1, we need x, a, k, sd_fun, boundary_sd, and intercept_sd
    if(is.null(sd_fun) || is.null(boundary_sd)){
      stop("sd_fun and boundary_sd are required for v1")
    }
    return(simulate_data_v1(x = x, a = a, k = k, sd_fun = sd_fun, boundary_sd = boundary_sd, intercept_sd = intercept_sd))

  } else if(type == "v2"){
    # For v2, we need x, a, sd_quasi, sd_exact, and intercept_sd
    return(simulate_data_v2(x = x, a = a, sd_quasi = sd_quasi, sd_exact = sd_exact, intercept_sd = intercept_sd))

  } else if(type == "np"){
    # For np, we need x, sd_np, and intercept_sd
    n <- length(x)
    # f <- rnorm(1, 0, sd_np) * x^2 + rnorm(1, 0, sd_np) * x + rnorm(1, 0, intercept_sd)

    # simulate some bspline basis
    f <- rnorm(1, 0, slope_sd)*x + rnorm(1, 0, sd_np) * bs(x, degree = 3, knots = quantile(x, probs = seq(0, 1, length.out = 5))) %*% rnorm(8)

    # for comparison, center f around some value
    f <- f - mean(f) + 3

    y <- rpois(lambda = exp(f), n = n)
    return(list(x = x, y = y, f = f))

  }
  else if(type == "iwp"){
    n <- length(x)

    x_min <- min(x); x_max <- max(x)
    knots <- seq(x_min, x_max, length.out = 50)
    pred_step = 1
    p = 2
    x_new_design <- BayesGP:::global_poly_helper(x = x, p = p)

    # Generate random weights for the linear functions
    beta_vec <- rnorm(n = p, mean = 0, sd = boundary_sd)

    if(sd_fun != 0){
      sd_fun <- sd_fun/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2)))

      prec_mat <- (1/sd_fun^2) * BayesGP:::compute_weights_precision_helper(knots)
      weights <- as.vector(LaplacesDemon::rmvnp(n = 1, mu = rep(0, ncol(prec_mat)), Omega = prec_mat))

      spline_new <- BayesGP:::local_poly_helper(knots = knots, refined_x = x, p = p)
      f <- x_new_design %*% beta_vec + as.vector(spline_new %*% weights)
    } else{
      f <- x_new_design %*% beta_vec
    }

    # for comparison, center f around some value
    f <- f - mean(f) + 3

    y <- rpois(lambda = exp(f), n = n)
    return(list(x = x, y = y, f = f))
  }
  else {
    stop("type not recognized")
  }
}
