library(boot)

grid_vec <- seq(0,1,by = 0.1)
seed <- 5; penalty <- 1;
datasets <- sim_dataset(
  pi = c(0.7, 0.15, 0.15),
  N = 300,
  sigma = 1,
  snr = NULL,
  seed = seed,
  pred_step = pred_step_selected,
  plot = T
)
L_matrix <- getLmat(datasets, grid_vec, order = 1, num_cores = 4, pred_step = pred_step_selected)
# alpha <- 0
# L_matrix_coarse <- (alpha/(alpha+nrow(L_matrix))) * L_matrix
# compute_estimates(L_matrix_coarse, penalty = 2)
compute_estimates(L_matrix, penalty = penalty)

# plot quantile of beta(a,a+k)
num_cluster <- length(grid_vec) - 1
alpha_vec <- seq(1,20, by = 1)
pbeta(q = 0.5, shape1 = alpha_vec, shape2 = num_cluster)

boot_mle <- boot(data = L_matrix,
                 stype = "i",
     statistic = function(data, indices) {
       data = data[indices, ]
       return(compute_estimates(data, penalty = penalty)$pi0_mle)
     },
     R = 3000
)

quantile(boot_mle$t)
sd(boot_mle$t)
quantile(boot_mle$t, 0.99)


if_mle <- empinf(boot.out = boot_mle, stype = "i", type = "reg")


boot_mom <- boot(data = L_matrix,
                 stype = "i",
     statistic = function(data, indices) {
       data = data[indices, ]
       return(compute_estimates(data)$pi0_mom)
     },
     R = 3000
)
if_mom <- empinf(boot.out = boot_mom, stype = "i", type = "reg")
plot(if_mle, if_mom, xlab = "MLE", ylab = "MOM", main = "Influence Functions")


## take a look at the influence functions for each observation
plot(if_mle ~ c(1:nrow(L_matrix)), type = "l",
     xlab = "Observation", ylab = "Influence Function", main = "")
points(if_mom ~ c(1:nrow(L_matrix)), type = "l", col = "red")

matplot(t(L_matrix[1:10,]), type = 'l')
matplot(t(L_matrix[288:298,]), type = 'l')


# pi0_grid <- seq(0.01, 0.99, by = 0.01)
# pi1_grid <- 1 - pi0_grid
# like_pi <- exp(L_matrix) %*% rbind(pi0_grid, pi1_grid)
# like_pi_vec <- colSums(log(like_pi))
# plot(like_pi_vec ~ pi0_grid)

like_diff <- apply(L_matrix, 1, function(x) {
  return(x[1] - max(x[-1]))
})

plot(like_diff ~ c(1:nrow(L_matrix)), type = "l",
     xlab = "Observation",
     ylab = "Difference in Likelihood")

selected_indices <- which(abs(like_diff) < 1)
compute_estimates(L_matrix[selected_indices,], penalty = 1)


## plot the change of est with respect to the subset of data
sorted_indices <- sort((like_diff), index.return = T)$ix
L_sort <- L_matrix[sorted_indices, ]
original_labels <- c(rep(1, 0.85*300), rep(0,0.15*300))
sort_labels <- original_labels[sorted_indices]
est_mom <- compute_estimates(L_matrix)$pi0_mom
est_mle <- compute_estimates(L_matrix, penalty = 1)$pi0_mle
true_proportion <- sum(original_labels) / 300
for (i in 1:298) {
  L_sort_subset <- L_sort[-(1:i), ]
  est <- compute_estimates(L_sort_subset, penalty = 1)
  est_mom <- c(est_mom, est$pi0_mom)
  est_mle <- c(est_mle, est$pi0_mle)
  true_proportion <- c(true_proportion, sum(sort_labels[-(1:i)]) / (300 - i) )
}

plot(est_mle ~ c(1:299), type = "o",
     ylim = c(0, 1),
     xlab = "Number of Observations Removed",
     ylab = "Estimate",
     main = "Estimate with respect to the subset of data")
points(est_mom ~ c(1:299), type = "o", col = "red")
points(true_proportion ~ c(1:299), type = "o", col = "blue")
legend("bottomright", legend = c("MLE", "MOM", "True Proportion"),
       col = c("black", "red", "blue"), lty = 1:1, cex = 0.8)

# take a look at the influential observations
mean(original_labels)
mean(original_labels[abs(if_mle) > quantile(abs(if_mle), 0.99)])

mean(original_labels)
mean(original_labels[abs(like_diff) > quantile(abs(like_diff), 0.99)])


