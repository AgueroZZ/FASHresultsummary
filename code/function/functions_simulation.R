### Write a function, to simulate a dataset with size n, given a function g
simulate_data <- function(n = 300, g, sd=0.1){
  x <- sort(runif(n, 0, 5))
  y <- g(x) + rnorm(n, sd=sd)
  return(data.frame(x=x, y=y,truef = g(x)))
}

### Write a function, to simulate a random function as g using cubic B-spline basis representation with random basis weights:
simulate_random_function <- function(n_basis = 5, sd_function = 1) {
  if(n_basis < 3) stop("n_basis must be greater than 3")
  # Define the range and knots for the B-spline basis
  x_min <- 0
  x_max <- 5
  
  # Generate equally spaced knots within the range [x_min, x_max]
  knots <- seq(x_min, x_max, length.out = n_basis - 3) 
  
  # Generate random weights for the basis functions
  weights <- rnorm((n_basis), mean = 0, sd = sd_function) # Adjusted to match the number of knots for cubic spline
  
  # Return a function that evaluates the spline at new x values
  function(x_new) {
    # Create the B-spline basis for the new x values using the predefined knots
    bspline_new <- bs(x_new, knots = knots, degree = 3, Boundary.knots = c(x_min, x_max))
    return(as.vector(bspline_new %*% weights))
  }
}

### simulate process: first draw a random function g, then draw a dataset from g
simulate_process <- function(n = 300, n_basis = 5, sd_fun = 1, sd = 0.1){
  g <- simulate_random_function(n_basis = n_basis, sd = sd_fun)
  return(simulate_data(n = n, g = g, sd = sd))
}

# plot(simulate_process(n = 300, n_basis = 3, sd_fun = 1, sd = 0.1))

