perform_rq <- function(data, tau_vec = seq(0.1, 0.9, by = 0.1), bootstrap = FALSE, bootstrap_B = 10000){
  beta_vec <- c()
  se_vec <- c()
  n <- nrow(data)

  if(!bootstrap){
    # Divide the data into non-overlapping blocks for each tau
    block_size <- floor(n / length(tau_vec))
    data_list <- split(data, rep(1:length(tau_vec), each = block_size, length.out = n))
  }
  else{
    data_list <- list()
    for(i in 1:length(tau_vec)){
      data_list[[i]] <- data[sample(1:n, n, replace = TRUE), ]
    }
  }

  # Perform quantile regression for each tau
  for (i in seq_along(tau_vec)) {
    tau <- tau_vec[i]
    mod <- rq(Y ~ G, data = data_list[[i]], tau = tau)
    summ_mod <- summary(mod, se = "boot")$coefficients
    beta_vec <- c(beta_vec, summ_mod[2, 1])
    se_vec <- c(se_vec, summ_mod[2, 2])
  }

  return(data.frame(tau = tau_vec, beta = beta_vec, se = se_vec))
}
