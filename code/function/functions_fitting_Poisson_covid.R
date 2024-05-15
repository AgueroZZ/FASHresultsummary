#############################################################
#############################################################
#############################################################
### Compute the log-marginal likelihood of the IWP model with O-Spline
compute_log_likelihood_ospline <- function(dataset, p, num_knots = 100, psd_iwp, pred_step, betaprec = 0.001){
  intercept <- -3 # assume this is the intercept
  x <- dataset$x
  y <- dataset$y
  weekdays <- dataset$weekdays
  knots <- unique(c(0,seq(min(x), max(x), length=num_knots)))
  # X1 <- BayesGP:::global_poly_helper(x, p = p)[,-1]
  X2 <- model.matrix(~ weekdays, data = dataset, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
  # X <- cbind(X1, X2)
  X <- X2
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
    intercept = intercept
  )

  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )

  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "Poisson_covid",
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
      intercept = intercept
    )

    tmbparams <- list(
      W = c(rep(0, (ncol(X))))
    )

    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = "Poisson_just_fixed_covid",
      silent = TRUE
    )
  }
  -ff$fn()
}

### Fit the IWP model with O-Spline
fit_ospline <- function(dataset, p, num_knots = 100, psd_iwp, pred_step, betaprec = 0.001){
  intercept <- -3 # assume this is the intercept
  condition_fitting <- 0 # 1 means the fitting is good, 1 means the fitting requires numerical adjustment
  x <- dataset$x
  y <- dataset$y
  weekdays <- dataset$weekdays
  knots <- unique(c(0,seq(min(x), max(x), length=num_knots)))
  X <- model.matrix(~ weekdays, data = dataset, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
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
    intercept = intercept,
    sigmaIWP = psd_iwp/sqrt((pred_step^((2 * p) - 1)) / (((2 * p) - 1) * (factorial(p - 1)^2)))
  )

  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )

  ff2 <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    DLL = "Poisson_covid",
    silent = TRUE
  )
  ff2$he <- function(w) numDeriv::jacobian(ff2$gr, w)

  opt <- nlminb(start = ff2$par, objective = ff2$fn, gradient = ff2$gr, hessian = ff2$he,
                control = list(eval.max = 20000, iter.max = 20000))

  prec_matrix <- Matrix::forceSymmetric(ff2$he(opt$par))
  smallest_eigen <- min(eigen(prec_matrix)$values)
  if(smallest_eigen < 0){
    condition_fitting <- 1
    ## Add a small value to the diagonal to make it positive definite
    prec_matrix <- prec_matrix + diag((abs(smallest_eigen) + 1e-5), nrow = ncol(prec_matrix))
  }

  mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)

  samps_coef <- LaplacesDemon::rmvnp(n = 8000, mu = mod$mean, Omega = as.matrix(mod$prec))
  samps_fitted <- as.matrix(B) %*% t(samps_coef[,1:ncol(B)])
  }
  else{
    tmbdat <- list(
      X = as(X,'dgTMatrix'),
      # Response
      y = y,
      # other known quantities
      betaprec = betaprec,
      intercept = intercept
    )

    tmbparams <- list(
      W = c(rep(0, (ncol(X))))
    )

    ff2 <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      DLL = "Poisson_just_fixed_covid",
      silent = TRUE
    )
    ff2$he <- function(w) numDeriv::jacobian(ff2$gr, w)

    opt <- nlminb(start = ff2$par, objective = ff2$fn, gradient = ff2$gr, hessian = ff2$he,
                  control = list(eval.max = 20000, iter.max = 20000))

    prec_matrix <- Matrix::forceSymmetric(ff2$he(opt$par))
    smallest_eigen <- min(eigen(prec_matrix)$values)
    if(smallest_eigen < 0){
      condition_fitting <- 1
      ## Add a small value to the diagonal to make it positive definite
      prec_matrix <- prec_matrix + diag((abs(smallest_eigen) + 0.0001), nrow = ncol(prec_matrix))
    }
    mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)

    samps_coef <- LaplacesDemon::rmvnp(n = 8000, mu = mod$mean, Omega = as.matrix(mod$prec))
    # Ensure output is a matrix, especially when dimension is 1
    if (!is.matrix(samps_coef)) {
      samps_coef <- matrix(samps_coef, ncol = length(samps_coef), nrow = 1)
    }
    # samps_fitted <- as.matrix(X1) %*% t(samps_coef)[1:ncol(X1),]
    samps_fitted <- NULL
  }
  log_likelihood <- -ff2$fn(opt$par) - 0.5 * determinant(as.matrix(ff2$he(opt$par)), logarithm = TRUE)$modulus + 0.5 * length(opt$par) * log(2 * pi)
  samps_fitted <- samps_fitted + intercept
  list(samps_coef = samps_coef, samps_fitted = samps_fitted, mod = mod, log_likelihood = log_likelihood, condition_fitting = condition_fitting,
       knots = knots, p = p)
}

### visualize the fitted result
visualize_fit <- function(x, fit_result, y = NULL, plot_samps = FALSE, original = FALSE){
  samps_fitted <- fit_result$samps_fitted
  if(original){
    samps_fitted <- exp(samps_fitted)
  }
  ## Produce the summary
  fitted_median <- apply(samps_fitted, 1, median)
  fitted_upper <- apply(samps_fitted, 1, quantile, probs = 0.975)
  fitted_lower <- apply(samps_fitted, 1, quantile, probs = 0.025)
  ## plot the fitted result
  if(is.null(y)){
    ## only plot the summary
    plot(x, fitted_median, col = 'red', type = "l")
    lines(x, fitted_upper, col = 'blue')
    lines(x, fitted_lower, col = 'blue')
    if(plot_samps){
      for(i in 1:ncol(samps_fitted)){
        lines(x, samps_fitted[,i], col = 'grey', lty = 2)
      }
    }
  }
  else{
    ## also plot the original data
    plot(x, y)
    lines(x, fitted_median, col = 'red')
    lines(x, fitted_upper, col = 'blue')
    lines(x, fitted_lower, col = 'blue')
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
compute_log_likelihood_ospline_seq <- function(dataset, p, num_knots, psd_iwp_vector, pred_step, betaprec = 0.001) {
  log_likelihoods <- sapply(psd_iwp_vector, function(psd_iwp) {
    compute_log_likelihood_ospline(
      dataset = dataset,
      p = p,
      num_knots = num_knots,
      psd_iwp = psd_iwp,
      pred_step = pred_step,
      betaprec = betaprec
    )
  })
  return(log_likelihoods)
}

### Fit a sequence of k models with ospline using different psd_iwp and different p and record the k marginal likelihoods
compute_log_likelihood_ospline_seq2 <- function(dataset, p_vec, num_knots, psd_iwp_vector, pred_step, betaprec = 0.001) {
  # Initialize a list to store the results for each p
  all_log_likelihoods <- c()

  # Loop over each value of p in p_vec
  for (p in p_vec) {
    # Compute log likelihoods for the current value of p
    log_likelihoods <- compute_log_likelihood_ospline_seq(
      dataset = dataset,
      p = p,
      num_knots = num_knots,
      psd_iwp_vector = psd_iwp_vector,
      pred_step = pred_step,
      betaprec = betaprec
    )
    all_log_likelihoods <- c(all_log_likelihoods, log_likelihoods)
  }

  return(all_log_likelihoods)
}


### Fit a sequence of k models based on a prior matrix: (when psd_iwp is varying)
fit_ospline_with_prior <- function(dataset, p, num_knots, prior_weight, pred_step, betaprec = 0.001, num_cores = 1) {
  # Filter out rows where prior weight is zero
  active_prior <- prior_weight[prior_weight$prior_weight > 0, ]

  # List to store results of fit_ospline for each active prior weight
  results_list <- vector("list", length = nrow(active_prior))

  # Determine whether to use parallel processing
  if (num_cores > 1 && .Platform$OS.type != "windows") {
    # Use mclapply for parallel processing on Unix-like systems
    results_list <- mclapply(seq_len(nrow(active_prior)), function(i) {
      psd_iwp <- active_prior$psd_iwp[i]
      fit_ospline(dataset, p, num_knots, psd_iwp, pred_step, betaprec)
    }, mc.cores = num_cores)
  }
  else {
    # Process sequentially
    for (i in seq_len(nrow(active_prior))) {
      psd_iwp <- active_prior$psd_iwp[i]
      results_list[[i]] <- fit_ospline(dataset, p, num_knots, psd_iwp, pred_step, betaprec)
    }
  }

  # Extract log likelihoods from results
  log_likelihoods <- sapply(results_list, function(result) result$log_likelihood)

  # Extract the fitting condition from results
  condition_fitting <- sapply(results_list, function(result) result$condition_fitting)

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
  list(fitted_results = results_list, posterior_weights = posterior_weight_matrix, condition_fitting = condition_fitting)
}

### Fit a sequence of k models based on a prior matrix: (when both psd_iwp and p are varying)
fit_ospline_with_prior2 <- function(dataset, num_knots, prior_weight, pred_step, betaprec = 0.001, num_cores = 1) {
  # Filter out rows where prior weight is zero
  active_prior <- prior_weight[prior_weight$prior_weight > 0, ]

  # Initialize lists to store the fitting results and data for posterior weights calculation
  results_list <- vector("list", length = nrow(active_prior))
  all_log_likelihoods <- numeric(length = nrow(active_prior))
  condition_fitting <- numeric(length = nrow(active_prior))
  # Loop through each entry in active_prior
  if (num_cores > 1 && .Platform$OS.type != "windows") {
    # Parallel processing for non-Windows platforms
    results_list <- mclapply(seq_len(nrow(active_prior)), function(i) {
      row <- active_prior[i, ]
      fit_result <- fit_ospline(dataset, row$p, num_knots, row$psd_iwp, pred_step, betaprec)
      fit_result
    }, mc.cores = num_cores)
    # Extract log likelihoods from results
    all_log_likelihoods <- sapply(results_list, function(result) result$log_likelihood)
    # Extract the fitting condition from results
    condition_fitting <- sapply(results_list, function(result) result$condition_fitting)

  } else {
    # Sequential processing
    for (i in seq_len(nrow(active_prior))) {
      print(rownames(active_prior)[i])
      row <- active_prior[i, ]
      fit_result <- fit_ospline(dataset, row$p, num_knots, row$psd_iwp, pred_step, betaprec)
      results_list[[i]] <- fit_result
      all_log_likelihoods[i] <- fit_result$log_likelihood
      condition_fitting[i] <- fit_result$condition_fitting
    }
  }

  # Calculate and normalize posterior weights
  max_log_likelihood <- max(all_log_likelihoods)
  stable_exp <- exp(all_log_likelihoods - max_log_likelihood)
  posterior_weights <- stable_exp * active_prior$prior_weight
  posterior_weights <- posterior_weights / sum(posterior_weights)

  # Create a data frame of posterior weights along with p and psd_iwp
  all_posterior_weights <- data.frame(p = active_prior$p, psd_iwp = active_prior$psd_iwp, posterior_weight = posterior_weights)

  return(list(fitted_results = results_list, posterior_weights = all_posterior_weights, condition_fitting = condition_fitting))
}


### aggregate result from fit_ospline_with_prior
aggregate_fit_with_prior <- function(x, fit_results_with_prior, original = FALSE) {
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

  if(original){
    aggregated_samples <- exp(aggregated_samples)
  }

  # Calculate mean and standard deviation across the sampled fits
  fitted_median <- apply(aggregated_samples, 1, median)
  fitted_upper <- apply(aggregated_samples, 1, quantile, probs = 0.975)
  fitted_lower <- apply(aggregated_samples, 1, quantile, probs = 0.025)
  summary_df <- data.frame(x = x, median = fitted_median, upper = fitted_upper, lower = fitted_lower)

  list(summary_df = summary_df, aggregated_samples = aggregated_samples)
}

### visualize result from fit_ospline_with_prior
visualize_fit_with_prior <- function(x, fit_results_with_prior, y = NULL, plot_samps = FALSE, original = FALSE) {
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

  if(original){
    aggregated_samples <- exp(aggregated_samples)
  }

  # Calculate mean and standard deviation across the sampled fits
  fitted_median <- apply(aggregated_samples, 1, median)
  fitted_upper <- apply(aggregated_samples, 1, quantile, probs = 0.975)
  fitted_lower <- apply(aggregated_samples, 1, quantile, probs = 0.025)

  # Create the plot
  plot_df <- data.frame(x = x, median = fitted_median, upper = fitted_upper, lower = fitted_lower)

  gg <- ggplot(plot_df, aes(x = x)) +
    geom_line(aes(y = median), color = 'red') +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = 'blue')

  if (!is.null(y)) {
    gg <- gg + geom_point(aes(y = y), size = 0.1)
  }

  if (plot_samps) {
    for (i in 1:ncol(aggregated_samples)) {
      gg <- gg + geom_line(aes(y = aggregated_samples[, i]), color = 'grey', linetype = "dotted") + theme_classic()
    }
  }

  print(gg)
}






