library(fashr)
library(tidyverse)
alpha <- 0.05
num_cores <- 32

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
load("./results/fash_fit2_all.RData")
datasets <- fash_fit2$fash_data$data_list
for (i in 1:length(datasets)) {
  datasets[[i]]$SE <- fash_fit2$fash_data$S[[i]]
}
cat("Total number of datasets: ", length(datasets), "\n")
fash_fit2$prior_weight
fash_fit2_update <- BF_update(fash_fit2, plot = FALSE)
fash_fit2_update$prior_weight


###############################################################
###############################################################
###############################################################
###############################################################
##################### 1. Compute min_lfsr_summary #############
###############################################################
###############################################################
###############################################################
smooth_var_refined = seq(0,15, by = 0.1)
min_lfsr_summary2 <- min_lfsr_summary(fash_fit2, num_cores = num_cores, smooth_var = smooth_var_refined)
save(min_lfsr_summary2, file = "./results/min_lfsr_summary2.RData")
min_lfsr_summary2_update <- min_lfsr_summary(fash_fit2_update, num_cores = num_cores, smooth_var = smooth_var_refined)
save(min_lfsr_summary2_update, file = "./results/min_lfsr_summary2_update.RData")

min_lfsr_summary2_deriv <- min_lfsr_summary(fash_fit2, num_cores = num_cores, smooth_var = smooth_var_refined, deriv = 1)
save(min_lfsr_summary2_deriv, file = "./results/min_lfsr_summary2_deriv.RData")
min_lfsr_summary2_update_deriv <- min_lfsr_summary(fash_fit2_update, num_cores = num_cores, smooth_var = smooth_var_refined, deriv = 1)
save(min_lfsr_summary2_update_deriv, file = "./results/min_lfsr_summary2_update_deriv.RData")
