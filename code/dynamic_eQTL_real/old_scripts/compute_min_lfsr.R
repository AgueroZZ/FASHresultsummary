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
datasets <- fash_fit1$fash_data$data_list
for (i in 1:length(datasets)) {
  datasets[[i]]$SE <- fash_fit1$fash_data$S[[i]]
}
cat("Total number of datasets: ", length(datasets), "\n")
original_weights <- fash_fit1$prior_weight
original_weights



###############################################################
###############################################################
###############################################################
###############################################################
########### 0. Update the fash obj based on BF ################
###############################################################
###############################################################
###############################################################
fash_fit1 <- BF_update(fash_fit1, plot = FALSE)



###############################################################
###############################################################
###############################################################
###############################################################
##################### 1. Compute min_lfsr_summary #############
###############################################################
###############################################################
###############################################################
smooth_var_refined = seq(0,15, by = 0.1)
min_lfsr_summary1 <- min_lfsr_summary(fash_fit1, num_cores = num_cores, smooth_var = smooth_var_refined)
save(min_lfsr_summary1, file = "./results/min_lfsr_summary1.RData")
