### Write a function, to simulate a dataset with size n, given a function g
simulate_data <- function(g, sd=0.1, snr = NULL){
  x <- 1:16
  if(!is.null(snr)){
    signal_var <- var(g(x))
    noise_var <- signal_var/snr
    sd <- sqrt(noise_var)
  }
  y <- g(x) + rnorm(n = length(x), sd = sd)
  return(data.frame(x=x, y=y,truef = g(x), sd = sd))
}

### Write a function, to simulate a random function as g using cubic B-spline basis representation with random basis weights:
simulate_nonlinear_function <- function(n_basis = 5, sd_function = 1, sd_poly = 0.1) {
  if(n_basis < 3) stop("n_basis must be greater than 3")
  # Define the range and knots for the B-spline basis
  x_min <- 0
  x_max <- 16

  # Generate equally spaced knots within the range [x_min, x_max]
  knots <- seq(x_min, x_max, length.out = n_basis - 3)

  # Generate random weights for the basis functions
  weights <- rnorm((n_basis), mean = 0, sd = sd_function) # Adjusted to match the number of knots for cubic spline

  # Generate random weights for the linear functions
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  beta1 <- rnorm(1, mean = 0, sd = sd_poly)

  # Return a function that evaluates the spline at new x values
  function(x_new) {
    # Create the B-spline basis for the new x values using the predefined knots
    bspline_new <- bs(x_new, knots = knots, degree = 3, Boundary.knots = c(x_min, x_max))
    # Return the function
    return(beta0 + beta1 * x_new + as.vector(bspline_new %*% weights))
  }
}

### Write a function, to simulate a random function as g using constant and linear function
simulate_linear_function <- function(sd_poly = 1){
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  beta1 <- rnorm(1, mean = 0, sd = sd_poly)
  function(x_new) {
    return(beta0 + beta1 * x_new)
  }
}

### Write a function, to simulate a random function as g using constant, linear and quadratic function
simulate_quadratic_function <- function(sd_poly = 1){
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  beta1 <- rnorm(1, mean = 0, sd = sd_poly)
  beta2 <- rnorm(1, mean = 0, sd = sd_poly)
  function(x_new) {
    return(beta0 + beta1 * x_new + beta2 * x_new^2)
  }
}

### Write a function, to simulate a nondynamic function as g using constant function
simulate_nondynamic_function <- function(sd_poly = 1){
  beta0 <- rnorm(1, mean = 0, sd = sd_poly)
  function(x_new) {
    return(beta0)
  }
}

### simulate process: first draw a random function g, then draw a dataset from g
simulate_process <- function(n_basis = 5, sd_fun = 1, sd = 0.1, snr = NULL, sd_poly = 0.1, type = "nonlinear"){
  if(type == "linear"){
    g <- simulate_linear_function(sd_poly = sd_poly)

  }
  else if(type == "quadratic"){
    g <- simulate_quadratic_function(sd_poly = sd_poly)
  }
  else if(type == "nonlinear") {
    g <- simulate_nonlinear_function(n_basis = n_basis, sd_function = sd_fun, sd_poly = sd_poly)
  }
  else if(type == "nondynamic") {
    g <- simulate_nondynamic_function(sd_poly = sd_poly)
  }
  else {
    stop("type must be one of 'linear', 'nonlinear', 'nondynamic'")
  }
  return(simulate_data(g = g, sd = sd, snr = snr))
}


