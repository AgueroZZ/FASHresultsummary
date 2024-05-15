#############################################################
#############################################################
#############################################################
### Compute the log-marginal likelihood of the sGP model with sB-Spline (yes!)
compute_log_likelihood_sBspline <- function(x, y, p, num_knots = 30, psd, pred_step, betaprec = 0.001, region){
  X1 <- BayesGP:::global_poly_helper(x, p = 2)
  X2 <- cbind(cos(2*pi*x), sin(2*pi*x))
  X <- cbind(X1, X2)
  P <- sGPfit::Compute_Q_sB(a = (2*pi), k = num_knots, region = region, accuracy = 0.001)
  B <- sGPfit::Compute_B_sB(x = x, a = (2*pi), k = num_knots, region = region)
  if(psd != 0){
  tmbdat <- list(
    X = as(X,'dgTMatrix'),
    B = as(B,'dgTMatrix'),
    P = as(P,'dgTMatrix'),
    logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
    # Response
    y = y,
    # other known quantities
    betaprec = betaprec,
    sigmaSmooth = psd/sGPfit:::compute_d_step_sGPsd(d = pred_step, a = (2*pi))
  )

  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )

  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "Poisson_causes",
    silent = TRUE
  )
  }
  else{
    tmbdat <- list(
      X = as(X,'dgTMatrix'),
      # Response
      y = y,
      # other known quantities
      betaprec = betaprec
    )

    tmbparams <- list(
      W = c(rep(0, (ncol(X))))
    )

    ff <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      random = "W",
      DLL = "Poisson_just_fixed_causes",
      silent = TRUE
    )
  }
  -ff$fn()
}

### Fit the sGP model with sB-Spline (yes!)
fit_sBspline <- function(x, y, num_knots = 30, psd, pred_step, betaprec = 0.001, region){
  X1 <- BayesGP:::global_poly_helper(x, p = 2)
  X2 <- cbind(cos(2*pi*x), sin(2*pi*x))
  X <- cbind(X1, X2)
  P <- sGPfit::Compute_Q_sB(a = (2*pi), k = num_knots, region = region, accuracy = 0.001)
  B <- sGPfit::Compute_B_sB(x = x, a = (2*pi), k = num_knots, region = region)
  if(psd != 0){
  tmbdat <- list(
    X = as(X,'dgTMatrix'),
    B = as(B,'dgTMatrix'),
    P = as(P,'dgTMatrix'),
    logPdet = as.numeric(determinant(P,logarithm = T)$modulus),
    # Response
    y = y,
    # other known quantities
    betaprec = betaprec,
    sigmaSmooth = psd/sGPfit:::compute_d_step_sGPsd(d = pred_step, a = (2*pi))
  )

  tmbparams <- list(
    W = c(rep(0, (ncol(X) + ncol(B))))
  )

  ff2 <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    DLL = "Poisson_causes",
    silent = TRUE
  )
  ff2$he <- function(w) numDeriv::jacobian(ff2$gr, w)

  opt <- nlminb(start = ff2$par, objective = ff2$fn, gradient = ff2$gr, hessian = ff2$he,
                control = list(eval.max = 20000, iter.max = 20000))

  prec_matrix <- Matrix::forceSymmetric(ff2$he(opt$par))
  mod = list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt)

  samps_coef <- LaplacesDemon::rmvnp(n = 8000, mu = mod$mean, Omega = as.matrix(mod$prec))
  samps_fitted <- as.matrix(B) %*% t(samps_coef[,1:ncol(B)]) + as.matrix(X) %*% t(samps_coef[,(ncol(B)+1):ncol(samps_coef)])
  }
  else{
    tmbdat <- list(
      X = as(X,'dgTMatrix'),
      # Response
      y = y,
      # other known quantities
      betaprec = betaprec
    )

    tmbparams <- list(
      W = c(rep(0, (ncol(X))))
    )

    ff2 <- TMB::MakeADFun(
      data = tmbdat,
      parameters = tmbparams,
      DLL = "Poisson_just_fixed_causes",
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

### visualize the fitted result (yes!)
visualize_fit <- function(x, fit_result, y = NULL, plot_samps = FALSE, original = FALSE){
  samps_fitted <- fit_result$samps_fitted
  if(original){
    samps_fitted <- exp(samps_fitted)
  }
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

### Fit a sequence of k models with sBspline using different psd and record the k marginal likelihoods (yes!)
compute_log_likelihood_sBspline_seq <- function(dataset, num_knots, psd_vector, pred_step, betaprec = 0.001, region){
  log_likelihoods <- sapply(psd_vector, function(psd) {
    compute_log_likelihood_sBspline(
      x = dataset$x,
      y = dataset$y,
      num_knots = num_knots,
      psd = psd,
      pred_step = pred_step,
      betaprec = betaprec,
      region = region
    )
  })
  return(log_likelihoods)
}

### Fit a sequence of k models based on a prior matrix: (when psd is varying)
fit_sBspline_with_prior <- function(x, y, num_knots, prior_weight, pred_step, betaprec = 0.001, region, num_cores = 1) {
  # Filter out rows where prior weight is zero
  active_prior <- prior_weight[prior_weight$prior_weight > 0, ]

  # List to store results of fit_ospline for each active prior weight
  results_list <- vector("list", length = nrow(active_prior))

  # Determine whether to use parallel processing
  if (num_cores > 1 && .Platform$OS.type != "windows") {
    # Use mclapply for parallel processing on Unix-like systems
    results_list <- mclapply(seq_len(nrow(active_prior)), function(i) {
      psd <- active_prior$psd[i]
      fit_sBspline(x, y, num_knots, psd, pred_step, betaprec, region)
     }, mc.cores = num_cores)
  }
  else {
    # Process sequentially
    for (i in seq_len(nrow(active_prior))) {
      psd <- active_prior$psd[i]
      results_list[[i]] <- fit_sBspline(x, y, num_knots, psd, pred_step, betaprec, region)
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
  posterior_weight_matrix <- cbind(active_prior$psd, posterior_weights)
  colnames(posterior_weight_matrix) <- c("psd", "posterior_weight")

  # Return both the list of fitted results and the matrix of posterior weights
  list(fitted_results = results_list, posterior_weights = posterior_weight_matrix)
}


### aggregate result from fit_ospline_with_prior
aggregate_fit_with_prior <- function(x, fit_results_with_prior, original = FALSE) {
  fit_results <- fit_results_with_prior$fitted_results
  posterior_weights <- fit_results_with_prior$posterior_weights[,"posterior_weight"]
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
  fitted_mean <- apply(aggregated_samples, 1, mean)
  fitted_sd <- apply(aggregated_samples, 1, sd)
  summary_df <- data.frame(x = x, mean = fitted_mean, upper = fitted_mean + 2 * fitted_sd, lower = fitted_mean - 2 * fitted_sd)

  list(summary_df = summary_df, aggregated_samples = aggregated_samples)
}

### visualize result from fit_ospline_with_prior
visualize_fit_with_prior <- function(x, fit_results_with_prior, y = NULL, plot_samps = FALSE, original = FALSE) {
  fit_results <- fit_results_with_prior$fitted_results
  posterior_weights <- fit_results_with_prior$posterior_weights[,"posterior_weight"]
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
  fitted_mean <- apply(aggregated_samples, 1, mean)
  fitted_sd <- apply(aggregated_samples, 1, sd)

  # Create the plot
  plot_df <- data.frame(x = x, mean = fitted_mean, upper = fitted_mean + 2 * fitted_sd, lower = fitted_mean - 2 * fitted_sd)

  gg <- ggplot(plot_df, aes(x = x)) +
    geom_line(aes(y = mean), color = 'red') +
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






