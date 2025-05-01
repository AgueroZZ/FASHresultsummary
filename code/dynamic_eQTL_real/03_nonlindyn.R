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
# load("./results/fash_fit1_all.RData")
load("./results/fash_fit2_all.RData")
datasets <- fash_fit2$fash_data$data_list
for (i in 1:length(datasets)) {
  datasets[[i]]$SE <- fash_fit2$fash_data$S[[i]]
}
cat("Total number of datasets: ", length(datasets), "\n")
original_weights <- fash_fit2$prior_weights


###############################################################
###############################################################
###############################################################
###############################################################
############# 2. Update the EB Result (Optional) ##############
###############################################################
###############################################################
###############################################################
# fash_fit1 <- BF_update(fash_fit1, plot = FALSE)
fash_fit2 <- BF_update(fash_fit2, plot = FALSE)


###############################################################
###############################################################
###############################################################
###############################################################
######## 9. Classification of nonlinear dynamics effect  ######
###############################################################
###############################################################
###############################################################

# Focus on the nonlinear dynamic eQTLs
alpha <- 0.05
test2 <- fdr_control(fash_fit2, alpha = alpha, plot = F)
selected_indices <- test2$fdr_results$index[test2$fdr_results$FDR <= alpha]

# Early (nonlin_dyn)
smooth_var_refined = seq(0,15, by = 0.1)
functional_early <- function(x){
  max(abs(x[smooth_var_refined <= 3])) - max(abs(x[smooth_var_refined > 3]))
}
testing_early_nonlin_dyn <- testing_functional(functional_early,
                                              lfsr_cal = function(x){mean(x <= 0)},
                                              fash = fash_fit2,
                                              indices = selected_indices,
                                              smooth_var = smooth_var_refined)
save(testing_early_nonlin_dyn, file = "./results/classify_nonlin_dyn_eQTLs_early.RData")

# Middle (nonlin_dyn)
functional_middle <- function(x){
  max(abs(x[smooth_var_refined <= 11 & smooth_var_refined >= 4])) - max(abs(x[smooth_var_refined > 11]), abs(x[smooth_var_refined < 4]))
}
testing_middle_nonlin_dyn <- testing_functional(functional_middle,
                                               lfsr_cal = function(x){mean(x <= 0)},
                                               fash = fash_fit2,
                                               indices = selected_indices,
                                               num_cores = num_cores,
                                               smooth_var = smooth_var_refined)
save(testing_middle_nonlin_dyn, file = "./results/classify_nonlin_dyn_eQTLs_middle.RData")

# Late (nonlin_dyn)
functional_late <- function(x){
  max(abs(x[smooth_var_refined >= 12])) - max(abs(x[smooth_var_refined < 12]))
}
testing_late_nonlin_dyn <- testing_functional(functional_late,
                                             lfsr_cal = function(x){mean(x <= 0)},
                                             fash = fash_fit2,
                                             indices = selected_indices,
                                             num_cores = num_cores,
                                             smooth_var = smooth_var_refined)
save(testing_late_nonlin_dyn, file = "./results/classify_nonlin_dyn_eQTLs_late.RData")

# Switch (nonlin_dyn)
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
testing_switch_nonlin_dyn <- testing_functional(functional_switch,
                                               lfsr_cal = function(x){mean(x <= 0)},
                                               fash = fash_fit2,
                                               indices = selected_indices,
                                               num_cores = num_cores,
                                               smooth_var = smooth_var_refined)
save(testing_switch_nonlin_dyn, file = "./results/classify_nonlin_dyn_eQTLs_switch.RData")



