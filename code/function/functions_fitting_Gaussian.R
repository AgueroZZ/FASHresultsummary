###########################################################
###########################################################
###########################################################
### For Precision matrix of the exact method: actual space
compute_Wp_cov <- function(x, p = 2){
  n <- length(x)
  result <- matrix(0, ncol = n, nrow = n)
  for (i in 1:n) {
    for (j in i:n) {
      ele <- 0:(p - 1)
      elements <- ((-1)^(p - 1 - ele)) * choose((2*p - 1),ele) * (x[i]^(2*p - 1 - ele)) * (x[j]^ele)
      term <- sum(elements)/factorial(2*p - 1)
      result[i,j] <- term 
    }
  }
  result <- Matrix::forceSymmetric(result)
  result
}
### Compute the log-marginal likelihood of the exact method with original-space
compute_log_likelihood_exact <- function(x, y, p, psd_iwp, pred_step, betaprec = 0.001, sd_gaussian = 0.1){
  C <- as.matrix(compute_Wp_cov(x, p = p))
  P = as(solve(C),'dgTMatrix')
  n <- length(x)
  X <- BayesGP:::global_poly_helper(x, p = p)
  B = Diagonal(n, x = 1)
  tmbdat <- list(
    X = as(X,'dgTMatrix'),
    # B = as(Diagonal(n, x = 1),'dgTMatrix'),
    B = B,
    P = P,
    # Response
    y = y,
    # other known quantities
    logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
    betaprec = betaprec,
    sigmaIWP = psd_iwp/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2))),
    sigmaGaussian = sd_gaussian
  )
  
  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )
  
  ff <- TMB::MakeADFun(
    random = "W",
    data = tmbdat,
    parameters = tmbparams,
    DLL = "Gaussian_theta_known",
    silent = TRUE
  )
  
  -ff$fn()
}
### Fit the exact method with original-space
fit_exact <- function(x, y, p, psd_iwp, pred_step, betaprec = 0.001, sd_gaussian = 0.1){
  C <- as.matrix(compute_Wp_cov(x, p = p))
  P = as(solve(C),'dgTMatrix')
  n <- length(x)
  X <- BayesGP:::global_poly_helper(x, p = p)
  B = Diagonal(n, x = 1)
  tmbdat <- list(
    X = as(X,'dgTMatrix'),
    # B = as(Diagonal(n, x = 1),'dgTMatrix'),
    B = B,
    P = P,
    # Response
    y = y,
    # other known quantities
    logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
    betaprec = betaprec,
    sigmaIWP = psd_iwp/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2))),
    sigmaGaussian = sd_gaussian
  )
  
  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )
  
  ff2 <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    DLL = "Gaussian_theta_known",
    silent = TRUE
  )
  ff2$he <- function(w) numDeriv::jacobian(ff2$gr, w)
  
  opt <- nlminb(start = ff2$par, objective = ff2$fn, gradient = ff2$gr, hessian = ff2$he, 
                control = list(eval.max = 20000, iter.max = 20000))
  
  prec_matrix <- Matrix::forceSymmetric(ff2$he(opt$par))
  mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)
  
  samps_coef <- LaplacesDemon::rmvnp(n = 8000, mu = mod$mean, Omega = as.matrix(mod$prec))
  samps_fitted <-as.matrix(B) %*% t(samps_coef[,1:ncol(B)]) + as.matrix(X) %*% t(samps_coef[,(ncol(B)+1):ncol(samps_coef)])
  
  log_likelihood <- -ff2$fn(opt$par) - 0.5 * determinant(as.matrix(ff2$he(opt$par)), logarithm = TRUE)$modulus + 0.5 * length(opt$par) * log(2 * pi)

  list(samps_coef = samps_coef, samps_fitted = samps_fitted, mod = mod, log_likelihood = log_likelihood)
  
}

#############################################################
#############################################################
#############################################################
### For Precision matrix of the exact method: augmented space
Compute_Ti <- function(x,p = 2,i){
  Ti <- matrix(0,nrow = p, ncol = p)
  delta <- diff(c(0,x))
  denom <- factorial(c(0:(p-1)))
  numr <- delta[i+1]^(0:(p-1))
  Ti[1,] <- numr/denom
  for (i in 2:p) {
    Ti[i,] <- c(rep(0,(i-1)),Ti[(i-1),((i-1):(p-1))])
  }
  Ti
}
Compute_Ci <- function(x, p = 2, i, is.cov = FALSE){
  delta <- diff(c(0,x))
  Result <- matrix(0,nrow = p, ncol = p)
  index <- i+1
  for (i in 1:p) {
    for (j in i:p) {
      Result[i,j] <- (delta[index]^(2*p + 1 - i - j))/((2*p + 1 - i - j)*factorial(p-i)*factorial(p-j))
    }
  }
  Result <- Matrix::forceSymmetric(Result)
  if(is.cov == T){
    return(Result)
  }
  else{
    round(solve(Result),digits = 5)
  }
}
Compute_Ai <- function(x, p = 2, i){
  Ci <- Compute_Ci(x,p,i)
  Ti <- Compute_Ti(x,p,i)
  Ai <- t(Ti) %*% Ci
  Ai <- Ai %*% Ti
  Ai
}
Compute_Bi <- function(x, p = 2, i){
  Ci <- Compute_Ci(x,p,i)
  Ti <- Compute_Ti(x,p,i)
  Bi <- -t(Ti) %*% Ci
  Bi
}
Compute_Aug_Wp_Prec <- function(x, p = 2) {
  n <- length(x)
  Blist <- list()
  AClist <- list()
  for (i in 1:(n - 1)) {
    AClist[[i]] <- Compute_Ai(x = x, i = i , p = p) + Compute_Ci(x = x, i = (i-1), p = p)
  }
  AClist[[n]] <- Compute_Ci(x = x, i = n - 1, p = p)
  for (i in 1:(n - 1)) {
    Blist[[i]] <- Compute_Bi(x = x, i = i, p = p)
  }
  Mlist <- list()
  M <- matrix(0, nrow = 0, ncol = n*p)
  for (i in 1:(n-1)) {
    Mlist[[i]] <- cbind(matrix(0,nrow = p, ncol = p * (i-1)), AClist[[i]], Blist[[i]], matrix(0,nrow = p, ncol = (p * (n-i-1))) )
    M <- rbind(M,Mlist[[i]])
  }
  M <- rbind(M,cbind(matrix(0,nrow = p, ncol = p * (n-1)), AClist[[n]]))
  M <- Matrix::forceSymmetric(M)
  as(as.matrix(M), "dgTMatrix")
}
Compute_design_Aug <- function(x, p = 2){
  Design <- Diagonal((p * length(x)), x = 0)
  diag(Design)[seq(1,nrow(Design), by = p)] <- 1
  as(as.matrix(Design[seq(1,nrow(Design), by = p),]), "dgTMatrix")
}
Compute_design_Aug_deriv <- function(x, p = 2, degree){
  Design <- Diagonal((p * length(x)), x = 0)
  diag(Design)[(seq(1,nrow(Design), by = p) + degree)] <- 1
  as(as.matrix(Design[(seq(1,nrow(Design), by = p) + degree),]), "dgTMatrix")
}
### Compute the log-marginal likelihood of the exact method with state-space
compute_log_likelihood_exact_aug <- function(x, y, p, psd_iwp, pred_step, betaprec = 0.001, sd_gaussian = 0.1){
  P = Compute_Aug_Wp_Prec(x, p = p)
  X = BayesGP:::global_poly_helper(x, p = p)
  B = Compute_design_Aug(x, p = p)
  tmbdat <- list(
    X = as(X,'dgTMatrix'),
    B = B,
    P = P,
    # Response
    y = y,
    # other known quantities
    logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
    betaprec = betaprec,
    sigmaIWP = psd_iwp/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2))),
    sigmaGaussian = sd_gaussian
  )
  
  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )
  
  ff <- TMB::MakeADFun(
    random = "W",
    data = tmbdat,
    parameters = tmbparams,
    DLL = "Gaussian_theta_known",
    silent = TRUE
  )
  
  -ff$fn()
}
### Fit the exact method with state-space
fit_exact_aug <- function(x, y, p, psd_iwp, pred_step, betaprec = 0.001, sd_gaussian = 0.1){
  P = Compute_Aug_Wp_Prec(x, p = p)
  X = BayesGP:::global_poly_helper(x, p = p)
  B = Compute_design_Aug(x, p = p)
  tmbdat <- list(
    X = as(X,'dgTMatrix'),
    B = B,
    P = P,
    # Response
    y = y,
    # other known quantities
    logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
    betaprec = betaprec,
    sigmaIWP = psd_iwp/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2))),
    sigmaGaussian = sd_gaussian
  )
  
  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )
  
  ff2 <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    DLL = "Gaussian_theta_known",
    silent = TRUE
  )
  ff2$he <- function(w) numDeriv::jacobian(ff2$gr, w)
  
  opt <- nlminb(start = ff2$par, objective = ff2$fn, gradient = ff2$gr, hessian = ff2$he, 
                control = list(eval.max = 20000, iter.max = 20000))
  
  prec_matrix <- Matrix::forceSymmetric(ff2$he(opt$par))
  mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)
  
  samps_coef <- LaplacesDemon::rmvnp(n = 8000, mu = mod$mean, Omega = as.matrix(mod$prec))
  samps_fitted <-as.matrix(B) %*% t(samps_coef[,1:ncol(B)]) + as.matrix(X) %*% t(samps_coef[,(ncol(B)+1):ncol(samps_coef)])
  
  log_likelihood <- -ff2$fn(opt$par) - 0.5 * determinant(as.matrix(ff2$he(opt$par)), logarithm = TRUE)$modulus + 0.5 * length(opt$par) * log(2 * pi)
  
  return(list(samps_coef = samps_coef, samps_fitted = samps_fitted, mod = mod, log_likelihood = log_likelihood))
}
  

#############################################################
#############################################################
#############################################################
### Compute the log-marginal likelihood of the IWP model with O-Spline
compute_log_likelihood_ospline <- function(x, y, p, num_knots = 100, psd_iwp, pred_step, betaprec = 0.001, sd_gaussian = 0.1){
  knots <- seq(min(x), max(x), length=num_knots)
  X <- BayesGP:::global_poly_helper(x, p = p)
  P <- BayesGP::compute_weights_precision_helper(knots)
  B <- BayesGP:::local_poly_helper(knots = knots, refined_x = x, p = p)
  if(psd_iwp != 0){
  tmbdat <- list(
    X = as(X,'dgTMatrix'),
    B = as(B,'dgTMatrix'),
    P = as(P,'dgTMatrix'),
    logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
    # Response
    y = y,
    # other known quantities
    betaprec = betaprec,
    sigmaIWP = psd_iwp/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2))),
    sigmaGaussian = sd_gaussian
  )
  
  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )
  
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "Gaussian_theta_known",
    silent = TRUE
  )
  }
  else{
    tmbdat <- list(
      X = as(X,'dgTMatrix'),
      # Response
      y = y,
      # other known quantities
      betaprec = betaprec,
      sigmaGaussian = sd_gaussian
    )
    
    tmbparams <- list(
      W = c(rep(0, (ncol(X))))
    )
    
    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = "Gaussian_just_fixed",
      silent = TRUE
    )
  }
  -ff$fn()
}

### Fit the IWP model with O-Spline
fit_ospline <- function(x, y, p, num_knots = 100, psd_iwp, pred_step, betaprec = 0.001, sd_gaussian = 0.1){
  knots <- seq(min(x), max(x), length=num_knots)
  X <- BayesGP:::global_poly_helper(x, p = p)
  P <- BayesGP::compute_weights_precision_helper(knots)
  B <- BayesGP:::local_poly_helper(knots = knots, refined_x = x, p = p)
  if(psd_iwp != 0){
  tmbdat <- list(
    X = as(X,'dgTMatrix'),
    B = as(B,'dgTMatrix'),
    P = as(P,'dgTMatrix'),
    logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
    # Response
    y = y,
    # other known quantities
    betaprec = betaprec,
    sigmaIWP = psd_iwp/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2))),
    sigmaGaussian = sd_gaussian
  )
  
  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )
  
  ff2 <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    DLL = "Gaussian_theta_known",
    silent = TRUE
  )
  ff2$he <- function(w) numDeriv::jacobian(ff2$gr, w)
  
  opt <- nlminb(start = ff2$par, objective = ff2$fn, gradient = ff2$gr, hessian = ff2$he, 
                control = list(eval.max = 20000, iter.max = 20000))
  
  prec_matrix <- Matrix::forceSymmetric(ff2$he(opt$par))
  mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)
  
  samps_coef <- LaplacesDemon::rmvnp(n = 8000, mu = mod$mean, Omega = as.matrix(mod$prec))
  samps_fitted <-as.matrix(B) %*% t(samps_coef[,1:ncol(B)]) + as.matrix(X) %*% t(samps_coef[,(ncol(B)+1):ncol(samps_coef)])
  }
  else{
    tmbdat <- list(
      X = as(X,'dgTMatrix'),
      # Response
      y = y,
      # other known quantities
      betaprec = betaprec,
      sigmaGaussian = sd_gaussian
    )
    
    tmbparams <- list(
      W = c(rep(0, (ncol(X))))
    )
    
    ff2 <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      DLL = "Gaussian_just_fixed",
      silent = TRUE
    )
    ff2$he <- function(w) numDeriv::jacobian(ff2$gr, w)
    
    opt <- nlminb(start = ff2$par, objective = ff2$fn, gradient = ff2$gr, hessian = ff2$he, 
                  control = list(eval.max = 20000, iter.max = 20000))
    
    prec_matrix <- Matrix::forceSymmetric(ff2$he(opt$par))
    mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)
    
    samps_coef <- LaplacesDemon::rmvnp(n = 8000, mu = mod$mean, Omega = as.matrix(mod$prec))
    # Ensure output is a matrix, especially when dimension is 1
    if (!is.matrix(samps_coef)) {
      samps_coef <- matrix(samps_coef, ncol = length(samps_coef), nrow = 1)
    }
    samps_fitted <- as.matrix(X) %*% t(samps_coef)
  }
  log_likelihood <- -ff2$fn(opt$par) - 0.5 * determinant(as.matrix(ff2$he(opt$par)), logarithm = TRUE)$modulus + 0.5 * length(opt$par) * log(2 * pi)

  list(samps_coef = samps_coef, samps_fitted = samps_fitted, mod = mod, log_likelihood = log_likelihood)
}

### visualize the fitted result
visualize_fit <- function(x, fit_result, y = NULL, plot_samps = FALSE){
  samps_fitted <- fit_result$samps_fitted
  ## Produce the summary
  fitted_mean <- apply(samps_fitted, 1, mean)
  fitted_sd <- apply(samps_fitted, 1, sd)
  ## plot the fitted result
  if(is.null(y)){
    ## only plot the summary
    plot(x, fitted_mean, col = 'red', type = "l")
    lines(x, (fitted_mean + (2*fitted_sd)), col = 'blue')
    lines(x, (fitted_mean - (2*fitted_sd)), col = 'blue')
    if(plot_samps){
      for(i in 1:ncol(samps_fitted)){
        lines(x, samps_fitted[,i], col = 'grey', lty = 2)
      }
    }
  }
  else{
    ## also plot the original data
    plot(x, y)
    lines(x, fitted_mean, col = 'red')
    lines(x, (fitted_mean + (2*fitted_sd)), col = 'blue')
    lines(x, (fitted_mean - (2*fitted_sd)), col = 'blue')
    if(plot_samps){
      for(i in 1:ncol(samps_fitted)){
        lines(x, samps_fitted[,i], col = 'grey', lty = 2)
      }
    }
  }
}


#################################################################
#################################################################
#################################################################
### Fit a sequence of k models with ospline using different psd_iwp and record the k marginal likelihoods
compute_log_likelihood_ospline_seq <- function(x, y, p, num_knots, psd_iwp_vector, pred_step, betaprec = 0.001, sd_gaussian = 0.1) {
  log_likelihoods <- sapply(psd_iwp_vector, function(psd_iwp) {
    compute_log_likelihood_ospline(
      x = x,
      y = y,
      p = p,
      num_knots = num_knots,
      psd_iwp = psd_iwp,
      pred_step = pred_step,
      betaprec = betaprec,
      sd_gaussian = sd_gaussian
    )
  })
  return(log_likelihoods)
}

### Fit a sequence of k models with ospline using different psd_iwp and different p and record the k marginal likelihoods
compute_log_likelihood_ospline_seq2 <- function(x, y, p_vec, num_knots, psd_iwp_vector, pred_step, betaprec = 0.001, sd_gaussian = 0.1) {
  # Initialize a list to store the results for each p
  all_log_likelihoods <- c()
  
  # Loop over each value of p in p_vec
  for (p in p_vec) {
    # Compute log likelihoods for the current value of p
    log_likelihoods <- compute_log_likelihood_ospline_seq(
      x = x,
      y = y,
      p = p,
      num_knots = num_knots,
      psd_iwp_vector = psd_iwp_vector,
      pred_step = pred_step,
      betaprec = betaprec,
      sd_gaussian = sd_gaussian
    )
    all_log_likelihoods <- c(all_log_likelihoods, log_likelihoods)
  }
  
  return(all_log_likelihoods)
}


### Fit a sequence of k models based on a prior matrix: (when psd_iwp is varying)
fit_ospline_with_prior <- function(x, y, p, num_knots, prior_weight, pred_step, betaprec = 0.001, sd_gaussian = 0.1, num_cores = 1) {
  # Filter out rows where prior weight is zero
  active_prior <- prior_weight[prior_weight$prior_weight > 0, ]
  
  # List to store results of fit_ospline for each active prior weight
  results_list <- vector("list", length = nrow(active_prior))
  
  # Determine whether to use parallel processing
  if (num_cores > 1 && .Platform$OS.type != "windows") {
    # Use mclapply for parallel processing on Unix-like systems
    results_list <- mclapply(seq_len(nrow(active_prior)), function(i) {
      psd_iwp <- active_prior$psd_iwp[i]
      fit_ospline(x, y, p, num_knots, psd_iwp, pred_step, betaprec, sd_gaussian)
    }, mc.cores = num_cores)
  }
  else {
    # Process sequentially
    for (i in seq_len(nrow(active_prior))) {
      psd_iwp <- active_prior$psd_iwp[i]
      results_list[[i]] <- fit_ospline(x, y, p, num_knots, psd_iwp, pred_step, betaprec, sd_gaussian)
    }
  }
  
  # Extract log likelihoods from results
  log_likelihoods <- sapply(results_list, function(result) result$log_likelihood)
  
  # Stabilize computation of posterior weights by centering log likelihoods
  max_log_likelihood <- max(log_likelihoods)
  stable_exp <- exp(log_likelihoods - max_log_likelihood)
  
  # Compute the posterior weights
  posterior_weights <- stable_exp * active_prior$prior_weight
  posterior_weights <- posterior_weights / sum(posterior_weights)
  
  # Create a matrix of posterior weights corresponding to active priors
  posterior_weight_matrix <- cbind(active_prior$psd_iwp, posterior_weights)
  colnames(posterior_weight_matrix) <- c("psd_iwp", "posterior_weight")
  
  # Return both the list of fitted results and the matrix of posterior weights
  list(fitted_results = results_list, posterior_weights = posterior_weight_matrix)
}

### Fit a sequence of k models based on a prior matrix: (when both psd_iwp and p are varying)
fit_ospline_with_prior2 <- function(x, y, num_knots, prior_weight, pred_step, betaprec = 0.001, sd_gaussian = 0.1, num_cores = 1) {
  # Filter out rows where prior weight is zero
  active_prior <- prior_weight[prior_weight$prior_weight > 0, ]
  
  # Initialize lists to store the fitting results and data for posterior weights calculation
  results_list <- vector("list", length = nrow(active_prior))
  all_log_likelihoods <- numeric(length = nrow(active_prior))
  
  # Loop through each entry in active_prior
  if (num_cores > 1 && .Platform$OS.type != "windows") {
    # Parallel processing for non-Windows platforms
    results_list <- mclapply(seq_len(nrow(active_prior)), function(i) {
      row <- active_prior[i, ]
      fit_result <- fit_ospline(x, y, row$p, num_knots, row$psd_iwp, pred_step, betaprec, sd_gaussian)
      fit_result
    }, mc.cores = num_cores)
    # Extract log likelihoods from results
    all_log_likelihoods <- sapply(results_list, function(result) result$log_likelihood)
  } else {
    # Sequential processing
    for (i in seq_len(nrow(active_prior))) {
      row <- active_prior[i, ]
      fit_result <- fit_ospline(x, y, row$p, num_knots, row$psd_iwp, pred_step, betaprec, sd_gaussian)
      results_list[[i]] <- fit_result
      all_log_likelihoods[i] <- fit_result$log_likelihood
    }
  }
  
  # Calculate and normalize posterior weights
  max_log_likelihood <- max(all_log_likelihoods)
  stable_exp <- exp(all_log_likelihoods - max_log_likelihood)
  posterior_weights <- stable_exp * active_prior$prior_weight
  posterior_weights <- posterior_weights / sum(posterior_weights)
  
  # Create a data frame of posterior weights along with p and psd_iwp
  all_posterior_weights <- data.frame(p = active_prior$p, psd_iwp = active_prior$psd_iwp, posterior_weight = posterior_weights)
  
  return(list(fitted_results = results_list, posterior_weights = all_posterior_weights))
}


### aggregate result from fit_ospline_with_prior
aggregate_fit_with_prior <- function(x, fit_results_with_prior, y = NULL, plot_samps = FALSE) {
  fit_results <- fit_results_with_prior$fitted_results
  posterior_weights <- fit_results_with_prior$posterior_weights[["posterior_weight"]]
  # Sample indices of the fit results based on posterior weights
  sampled_indices <- sample(seq_along(fit_results), size = 8000, replace = TRUE, prob = posterior_weights)
  
  # Tabulate the frequency of each index
  tabulated_indices <- table(sampled_indices)
  
  # Initialize an empty list to collect samples
  sampled_fits_list <- vector("list", length = length(tabulated_indices))
  
  # Retrieve and store the required number of columns from each sampled fit's samps_fitted
  names(tabulated_indices) <- as.integer(names(tabulated_indices))  # Ensure names are integer
  for (i in seq_along(tabulated_indices)) {
    idx <- as.integer(names(tabulated_indices)[i])
    count <- tabulated_indices[[i]]
    sampled_fits_list[[i]] <- fit_results[[idx]]$samps_fitted[, 1:count]
  }
  
  # Combine all samples into one matrix
  aggregated_samples <- do.call(cbind, sampled_fits_list)
  
  # Calculate mean and standard deviation across the sampled fits
  fitted_mean <- apply(aggregated_samples, 1, mean)
  fitted_sd <- apply(aggregated_samples, 1, sd)
  summary_df <- data.frame(x = x, mean = fitted_mean, upper = fitted_mean + 2 * fitted_sd, lower = fitted_mean - 2 * fitted_sd)

  list(summary_df = summary_df, aggregated_samples = aggregated_samples)
}

### visualize result from fit_ospline_with_prior
visualize_fit_with_prior <- function(x, fit_results_with_prior, y = NULL, plot_samps = FALSE, truef = NULL) {
  fit_results <- fit_results_with_prior$fitted_results
  posterior_weights <- fit_results_with_prior$posterior_weights[["posterior_weight"]]
  # Sample indices of the fit results based on posterior weights
  sampled_indices <- sample(seq_along(fit_results), size = 8000, replace = TRUE, prob = posterior_weights)
  
  # Tabulate the frequency of each index
  tabulated_indices <- table(sampled_indices)
  
  # Initialize an empty list to collect samples
  sampled_fits_list <- vector("list", length = length(tabulated_indices))
  
  # Retrieve and store the required number of columns from each sampled fit's samps_fitted
  names(tabulated_indices) <- as.integer(names(tabulated_indices))  # Ensure names are integer
  for (i in seq_along(tabulated_indices)) {
    idx <- as.integer(names(tabulated_indices)[i])
    count <- tabulated_indices[[i]]
    sampled_fits_list[[i]] <- fit_results[[idx]]$samps_fitted[, 1:count]
  }
  
  # Combine all samples into one matrix
  aggregated_samples <- do.call(cbind, sampled_fits_list)
  
  # Calculate mean and standard deviation across the sampled fits
  fitted_mean <- apply(aggregated_samples, 1, mean)
  fitted_sd <- apply(aggregated_samples, 1, sd)
  
  # Create the plot
  plot_df <- data.frame(x = x, mean = fitted_mean, upper = fitted_mean + 2 * fitted_sd, lower = fitted_mean - 2 * fitted_sd)
  
  gg <- ggplot(plot_df, aes(x = x)) +
    geom_line(aes(y = mean), color = 'red') +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = 'blue') +
    theme_classic()
  
  if (!is.null(y)) {
    gg <- gg + geom_point(aes(y = y), size = 0.1)
  }
  
  if(!is.null(truef)){
    gg <- gg + geom_line(aes(y = truef), color = 'black', linetype = "dashed")
  }
  
  if (plot_samps) {
    for (i in 1:ncol(aggregated_samples)) {
      gg <- gg + geom_line(aes(y = aggregated_samples[, i]), color = 'grey', linetype = "dotted") + theme_classic()
    }
  }
  
  print(gg)
}






