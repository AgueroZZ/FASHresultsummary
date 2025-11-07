compute_estimates <- function(L_matrix, penalty = 1){
  if (penalty > 1) {
    prior_null <- matrix(0, nrow = floor(penalty - 1), ncol = ncol(L_matrix))
    prior_null[, 1] <- 1  # Prior mass on the first grid point
    L_matrix_original <- rbind(exp(L_matrix), prior_null)
    fit.sqp <- mixsqp::mixsqp(
      L = L_matrix_original,
      log = FALSE,
      control = list(verbose = FALSE)
    )
    pi0_mle <- fit.sqp$x[1]

  } else{
    fit.sqp <- mixsqp::mixsqp(L = L_matrix, log = TRUE, control = list(verbose = FALSE))
    pi0_mle <- fit.sqp$x[1]
  }
  pi0_mom <- sum(apply(L_matrix, 1, function(x) which.max(x) == 1))/nrow(L_matrix)
  return(data.frame(pi0_mle = pi0_mle, pi0_mom = pi0_mom))
}

### Take a look at the first model
L_matrix <- fash_fit1$L_matrix

### The MLE estimate
penalty_select <- 1
pi_0_est <- compute_estimates(L_matrix, penalty = penalty_select)
pi_0_est
# pi0_mle   pi0_mom
# 1       0 0.4384977


boot_mle <- boot(data = L_matrix,
                 stype = "i",
                 statistic = function(data, indices, sub_samps = 15000) {
                   indices <- indices[sample(length(indices), sub_samps, replace = FALSE)]
                   data = data[indices, ]
                   return(compute_estimates(data, penalty = penalty_select)$pi0_mle)
                 },
                 R = 100,
                 parallel = "multicore",
                 ncpus = 2
)
quantile(boot_mle$t)


boot_mle2 <- boot(data = L_matrix,
                 stype = "i",
                 statistic = function(data, indices, sub_samps = 15000) {
                   indices <- indices[sample(length(indices), sub_samps, replace = FALSE)]
                   data = data[indices, ]
                   return(compute_estimates(data, penalty = 50)$pi0_mle)
                 },
                 R = 100,
                 parallel = "multicore",
                 ncpus = 2
)
quantile(boot_mle2$t)


num_cluster <- ncol(L_matrix) - 1
alpha_vec <- seq(1,100, by = 1)
qbeta_quant <- pbeta(q = 0.5, shape1 = alpha_vec, shape2 = num_cluster)
plot(qbeta_quant ~ alpha_vec, type = "l", xlab = "prob", ylab = "alpha", main = "Exceeding Prob")

