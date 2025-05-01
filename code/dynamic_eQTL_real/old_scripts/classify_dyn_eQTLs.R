library(fashr)
library(tidyverse)
alpha <- 0.05
num_cores <- 16
###############################################################
###############################################################
###############################################################
###############################################################
##################### 0. Loading ##############################
###############################################################
###############################################################
###############################################################
log_prec <- seq(0,10, by = 0.2)
fine_grid <- sort(c(0, exp(-0.5*log_prec)))
load("./results/fash_fit1_all.RData")
# load("./results/fash_fit2_all.RData")
datasets <- fash_fit1$fash_data$data_list
for (i in 1:length(datasets)) {
  datasets[[i]]$SE <- fash_fit1$fash_data$S[[i]]
}
cat("Total number of datasets: ", length(datasets), "\n")


###############################################################
###############################################################
###############################################################
###############################################################
############# 2. Update the EB Result (Optional) ##############
###############################################################
###############################################################
###############################################################
# Update the EB result:
fash_fit1 <- BF_update(fash_fit1, plot = FALSE)



###############################################################
###############################################################
###############################################################
###############################################################
######## 9. Classification of linear dynamics  ############
###############################################################
###############################################################
###############################################################
# Early (dyn)
smooth_var_refined = seq(0,15, by = 0.1)
functional_early <- function(x){
  max(abs(x[smooth_var_refined <= 3])) - max(abs(x[smooth_var_refined > 3]))
}
testing_early_dyn <- testing_functional(functional_early,
                                              lfsr_cal = function(x){mean(x <= 0)},
                                              fash = fash_fit1,
                                              indices = 1:length(datasets),
                                              smooth_var = smooth_var_refined)
save(testing_early_dyn, file = "./results/classify_dyn_eQTLs_early.RData")
# plot(testing_early_dyn$cfsr, type = "l", xlab = "Number of datasets", ylab = "FSR", main = "Early dyn")
# abline(h = alpha, col = "red", lty = 2)
# early_dyn_indices <- testing_early_dyn$indices[testing_early_dyn$cfsr <= alpha]
# par(mfrow = c(2,2))
# for (i in 1:4) {
#   selected_index <- sample(early_dyn_indices, 1)
#   fitted_result <- predict(fash_fit1,
#                            index = selected_index,
#                            smooth_var = seq(0, 15, by = 0.1))
#   plot(
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y,
#     pch = 20,
#     col = "black",
#     xlab = "Time",
#     ylab = "Effect Est",
#     main = paste0(names(datasets)[selected_index]),
#     ylim = c(
#       min(datasets[[selected_index]]$y - 2 * datasets[[selected_index]]$SE),
#       max(datasets[[selected_index]]$y + 2 * datasets[[selected_index]]$SE)
#     )
#   )
#   arrows(
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y - 2 * datasets[[selected_index]]$SE,
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y + 2 * datasets[[selected_index]]$SE,
#     length = 0.05,
#     angle = 90,
#     code = 3,
#     col = "black"
#   )
#   lines(fitted_result$x,
#         fitted_result$mean,
#         col = "red",
#         lwd = 2)
#   abline(h = 0, lty = 2, col = "blue")
#   polygon(
#     c(fitted_result$x, rev(fitted_result$x)),
#     c(fitted_result$lower, rev(fitted_result$upper)),
#     col = rgb(1, 0, 0, 0.3),
#     border = NA
#   )
# }
# par(mfrow = c(1,1))


# Middle (dyn)
functional_middle <- function(x){
  max(abs(x[smooth_var_refined <= 11 & smooth_var_refined >= 4])) - max(abs(x[smooth_var_refined > 11]), abs(x[smooth_var_refined < 4]))
}
testing_middle_dyn <- testing_functional(functional_middle,
                                               lfsr_cal = function(x){mean(x <= 0)},
                                               fash = fash_fit1,
                                               indices = 1:length(datasets),
                                               num_cores = num_cores,
                                               smooth_var = smooth_var_refined)
save(testing_middle_dyn, file = "./results/classify_dyn_eQTLs_middle.RData")
# plot(testing_middle_dyn$cfsr, type = "l", xlab = "Number of datasets", ylab = "FSR", main = "Middle dyn")
# abline(h = alpha, col = "red", lty = 2)
# middle_dyn_indices <- testing_middle_dyn$indices[testing_middle_dyn$cfsr <= alpha]
# par(mfrow = c(2,2))
# for (i in 1:4) {
#   selected_index <- sample(middle_dyn_indices, 1)
#   fitted_result <- predict(fash_fit1,
#                            index = selected_index,
#                            smooth_var = seq(0, 15, by = 0.1))
#   plot(
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y,
#     pch = 20,
#     col = "black",
#     xlab = "Time",
#     ylab = "Effect Est",
#     main = paste0(names(datasets)[selected_index]),
#     ylim = c(
#       min(datasets[[selected_index]]$y - 2 * datasets[[selected_index]]$SE),
#       max(datasets[[selected_index]]$y + 2 * datasets[[selected_index]]$SE)
#     )
#   )
#   arrows(
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y - 2 * datasets[[selected_index]]$SE,
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y + 2 * datasets[[selected_index]]$SE,
#     length = 0.05,
#     angle = 90,
#     code = 3,
#     col = "black"
#   )
#   lines(fitted_result$x,
#         fitted_result$mean,
#         col = "red",
#         lwd = 2)
#   abline(h = 0, lty = 2, col = "blue")
#   polygon(
#     c(fitted_result$x, rev(fitted_result$x)),
#     c(fitted_result$lower, rev(fitted_result$upper)),
#     col = rgb(1, 0, 0, 0.3),
#     border = NA
#   )
# }
# par(mfrow = c(1,1))

# Late (dyn)
functional_late <- function(x){
  max(abs(x[smooth_var_refined >= 12])) - max(abs(x[smooth_var_refined < 12]))
}
testing_late_dyn <- testing_functional(functional_late,
                                             lfsr_cal = function(x){mean(x <= 0)},
                                             fash = fash_fit1,
                                             indices = 1:length(datasets),
                                             num_cores = num_cores,
                                             smooth_var = smooth_var_refined)
save(testing_late_dyn, file = "./results/classify_dyn_eQTLs_late.RData")
# plot(testing_late_dyn$cfsr, type = "l", xlab = "Number of datasets", ylab = "FSR", main = "Late dyn")
# abline(h = alpha, col = "red", lty = 2)
# late_dyn_indices <- testing_late_dyn$indices[testing_late_dyn$cfsr <= alpha]
# par(mfrow = c(2,2))
# for (i in 1:4) {
#   selected_index <- sample(late_dyn_indices, 1)
#   fitted_result <- predict(fash_fit1,
#                            index = selected_index,
#                            smooth_var = seq(0, 15, by = 0.1))
#   plot(
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y,
#     pch = 20,
#     col = "black",
#     xlab = "Time",
#     ylab = "Effect Est",
#     main = paste0(names(datasets)[selected_index]),
#     ylim = c(
#       min(datasets[[selected_index]]$y - 2 * datasets[[selected_index]]$SE),
#       max(datasets[[selected_index]]$y + 2 * datasets[[selected_index]]$SE)
#     )
#   )
#   arrows(
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y - 2 * datasets[[selected_index]]$SE,
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y + 2 * datasets[[selected_index]]$SE,
#     length = 0.05,
#     angle = 90,
#     code = 3,
#     col = "black"
#   )
#   lines(fitted_result$x,
#         fitted_result$mean,
#         col = "red",
#         lwd = 2)
#   abline(h = 0, lty = 2, col = "blue")
#   polygon(
#     c(fitted_result$x, rev(fitted_result$x)),
#     c(fitted_result$lower, rev(fitted_result$upper)),
#     col = rgb(1, 0, 0, 0.3),
#     border = NA
#   )
# }
# par(mfrow = c(1,1))

# Switch (dyn)
switch_threshold <- 0.25
functional_switch <- function(x){
  # compute the radius of x, measured by deviation from 0 from below and from above
  x_pos <- x[x > 0]
  x_neg <- x[x < 0]
  if(length(x_pos) == 0 || length(x_neg) == 0){
    return(0)
  }
  min(max(abs(x_pos)), max(abs(x_neg))) - switch_threshold
}
testing_switch_dyn <- testing_functional(functional_switch,
                                               lfsr_cal = function(x){mean(x <= 0)},
                                               fash = fash_fit1,
                                               indices = 1:length(datasets),
                                               num_cores = num_cores,
                                               smooth_var = smooth_var_refined)
save(testing_switch_dyn, file = "./results/classify_dyn_eQTLs_switch.RData")
# plot(testing_switch_dyn$cfsr, type = "l", xlab = "Number of datasets", ylab = "FSR", main = "Switch dyn")
# abline(h = alpha, col = "red", lty = 2)
# switch_dyn_indices <- testing_switch_dyn$indices[testing_switch_dyn$cfsr <= alpha]
# par(mfrow = c(2,2))
# for (i in 1:4) {
#   selected_index <- sample(switch_dyn_indices, 1)
#   fitted_result <- predict(fash_fit1,
#                            index = selected_index,
#                            smooth_var = seq(0, 15, by = 0.1))
#   plot(
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y,
#     pch = 20,
#     col = "black",
#     xlab = "Time",
#     ylab = "Effect Est",
#     main = paste0(names(datasets)[selected_index]),
#     ylim = c(
#       min(datasets[[selected_index]]$y - 2 * datasets[[selected_index]]$SE),
#       max(datasets[[selected_index]]$y + 2 * datasets[[selected_index]]$SE)
#     )
#   )
#   arrows(
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y - 2 * datasets[[selected_index]]$SE,
#     datasets[[selected_index]]$x,
#     datasets[[selected_index]]$y + 2 * datasets[[selected_index]]$SE,
#     length = 0.05,
#     angle = 90,
#     code = 3,
#     col = "black"
#   )
#   lines(fitted_result$x,
#         fitted_result$mean,
#         col = "red",
#         lwd = 2)
#   abline(h = 0, lty = 2, col = "blue")
#   polygon(
#     c(fitted_result$x, rev(fitted_result$x)),
#     c(fitted_result$lower, rev(fitted_result$upper)),
#     col = rgb(1, 0, 0, 0.3),
#     border = NA
#   )
# }
# par(mfrow = c(1,1))


# Save the classification results
save(testing_early_dyn, testing_middle_dyn, testing_late_dyn, testing_switch_dyn, file = "./results/classify_dyn_eQTLs.RData")

