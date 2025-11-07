library(BayesGP)
library(TMB)
library(Matrix)
library(splines)
library(parallel)
library(ggplot2)
library(reshape2)
library(mixsqp)
library(tidyverse)

cpp_dir <- paste0(getwd(), "/cpp")
fig_dir <- paste0(dirname(getwd()), "/output/simulation")
result_dir <- paste0(dirname(getwd()), "/output/simulation/figure")
function_dir <- paste0(getwd(), "/function")
source(paste0(function_dir, "/functions_fitting.R"))
source(paste0(function_dir, "/functions_simulation.R"))
compile(paste0(cpp_dir, "/Gaussian_theta_known.cpp"))
compile(paste0(cpp_dir, "/Gaussian_just_fixed.cpp"))
dyn.load(TMB::dynlib(paste0(cpp_dir, "/Gaussian_theta_known")))
dyn.load(TMB::dynlib(paste0(cpp_dir, "/Gaussian_just_fixed")))

# Number of cores to use
num_cores <- 4

### Suppose we have 48 datasets in total:
## 16 datasets are simulated with n_basis = 5, sd_fun = 1 (category A)
## and 16 datasets are simulated with n_basis = 10, sd_fun = 1 (category B)
## and 16 datasets are simulated with n_basis = 20, sd_fun = 1 (category C)

# Generate indices of groups
group_indices <- rep(1:3, each = 16)

# Generate datasets in parallel
datasets <- mclapply(1:48, function(i) {
  set.seed(i)
  if (i <= 16) {
    n_basis <- 5
    sd_fun <- 1
  } else if (i <= 32) {
    n_basis <- 10
    sd_fun <- 1
  } else {
    n_basis <- 20
    sd_fun <- 1
  }
  simulate_process(n = 100, n_basis = n_basis, sd_fun = sd_fun, sd = 0.1)
}, mc.cores = num_cores)

# visualize them:
pdf(paste0(fig_dir, "/raw_datasets_catA.pdf"), height = 10, width = 10)
par(mfrow = c(4, 4), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (i in 1:16) {
  plot(datasets[[i]]$x, datasets[[i]]$y, type = 'p',  # Change to 'p' for points
       main = paste("Dataset", i), xlab = "x", ylab = "y",
       cex = 0.1,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  lines(datasets[[i]]$x, datasets[[i]]$truef, col = "red")
}
dev.off()
pdf(paste0(fig_dir, "/raw_datasets_catB.pdf"), height = 10, width = 10)
par(mfrow = c(4, 4), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (i in 17:32) {
  plot(datasets[[i]]$x, datasets[[i]]$y, type = 'p',  # Change to 'p' for points
       main = paste("Dataset", i), xlab = "x", ylab = "y",
       cex = 0.1,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  lines(datasets[[i]]$x, datasets[[i]]$truef, col = "red")
}
dev.off()
pdf(paste0(fig_dir, "/raw_datasets_catC.pdf"), height = 10, width = 10)
par(mfrow = c(4, 4), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (i in 33:48) {
  plot(datasets[[i]]$x, datasets[[i]]$y, type = 'p',  # Change to 'p' for points
       main = paste("Dataset", i), xlab = "x", ylab = "y",
       cex = 0.1,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  lines(datasets[[i]]$x, datasets[[i]]$truef, col = "red")
}
dev.off()
par(mfrow = c(1, 1))

# Compute L_vecs for each dataset in parallel
set.seed(123)
p_vec <- 2
# psd_iwp_vec <- seq(0,5, by = 0.1)
log_prec <- unique(sort(c(Inf, seq(-1,1, by = 0.1), seq(-5,-1, by = 0.5), seq(1,5, by = 0.5), seq(-10,-5, by = 1), seq(5,10, by = 1)), decreasing = T))
psd_iwp_vec <- 1/exp(.5*log_prec)
L_vecs <- mclapply(datasets, function(dataset) {
  compute_log_likelihood_ospline_seq2(
    x = dataset$x,
    y = dataset$y,
    p = p_vec,
    num_knots = 50,
    psd_iwp_vector = psd_iwp_vec,
    pred_step = 1,
    betaprec = 0.001,
    sd_gaussian = 0.1
  )
}, mc.cores = num_cores)
L_matrix <- do.call(rbind, L_vecs)
save(L_matrix, file = paste0(result_dir, "/L_matrix.rda"))
load(paste0(result_dir, "/L_matrix.rda"))
max_indices <- apply(L_matrix, 1, which.max)

# visualize the group_indices and the max_indices
plot(group_indices, max_indices, xlab = "Group indices", ylab = "Max indices")

### Let's learn the weights through EB:
fit.sqp <- mixsqp(L = L_matrix, log = TRUE)
numiter <- nrow(fit.sqp$progress)
plot(1:numiter,fit.sqp$progress$objective,type = "b",
     pch = 20,lwd = 2,xlab = "SQP iteration",
     ylab = "objective",xaxp = c(1,numiter,numiter - 1))
fit.sqp$x ## the optimal weights
prior_weight <- data.frame(p = rep(p_vec, each = length(psd_iwp_vec)), psd_iwp = psd_iwp_vec, prior_weight = fit.sqp$x)

### Visualize the fit:
fit_result_final <- fit_ospline_with_prior2(num_cores = 2, x = datasets[[1]]$x, y = datasets[[1]]$y, num_knots = 50, prior_weight = prior_weight, betaprec = 0.001, sd_gaussian = 0.1, pred_step = 1)
visualize_fit_with_prior(fit_results_with_prior = fit_result_final,
                         x = datasets[[1]]$x, y = datasets[[1]]$y, truef = datasets[[1]]$truef)


### Loop through the dataset and produce matrix posterior_weights
# Predefine the number of datasets
num_datasets <- length(datasets)
num_weights <- sum(prior_weight$prior_weight != 0)
posterior_weights_matrix <- matrix(nrow = num_datasets, ncol = num_weights)

# Loop through each dataset and perform fitting
fitted_datasets <- list()
for (i in seq_along(datasets)) {
  dataset <- datasets[[i]]
  fit_result_final <- fit_ospline_with_prior2(
    num_cores = num_cores,
    x = dataset$x,
    y = dataset$y,
    num_knots = 50,
    prior_weight = prior_weight,
    betaprec = 0.001,
    sd_gaussian = 0.1,
    pred_step = 1
  )
  fitted_datasets[[i]] <- aggregate_fit_with_prior(x = dataset$x, fit_results_with_prior = fit_result_final)$summary_df
  # Assume posterior weights are in the same order as in 'prior_weight$psd_iwp'
  posterior_weights_matrix[i, ] <- fit_result_final$posterior_weights[, "posterior_weight"]
}
colnames(posterior_weights_matrix) <- paste(as.character(fit_result_final$posterior_weights[, "p"]),
                                            as.character(fit_result_final$posterior_weights[, "psd_iwp"]), sep = "_")
save(posterior_weights_matrix, file = paste0(result_dir, "/posterior_weights_matrix.rda"))
save(fitted_datasets, file = paste0(result_dir, "/fitted_datasets.rda"))
load(paste0(result_dir, "/posterior_weights_matrix.rda"))
load(paste0(result_dir, "/fitted_datasets.rda"))

# Visualize the posterior weights
posterior_weights_df <- as.data.frame(posterior_weights_matrix)
posterior_weights_df$id <- 1:nrow(posterior_weights_df)
melted_data <- melt(posterior_weights_df, id.vars = "id")
melted_data$variable2 <- sub("_.*", "", melted_data$variable)
melted_data$variable3 <- as.factor(round(as.numeric(sub("*._", "", melted_data$variable)), 3))

ggplot(melted_data, aes(x = as.factor(id), y = value, fill = variable3)) +
  geom_bar(stat = "identity") +
  labs(x = "Observation ID", y = "Weight", fill = "PSD") +
  theme_minimal() +
  ggtitle("Structure Plot of Posterior Weights") +
  coord_flip() 
ggsave(paste0(fig_dir, "/posterior_weights_structure_plot.pdf"), height = 10, width = 10)

# Normalize the data (optional but recommended)
wide_data <- scale(posterior_weights_matrix)
# Perform hierarchical clustering
d <- dist(wide_data, method = "euclidean")  # Distance matrix
fit <- hclust(d, method = "ward.D2")  # Clustering
clusters <- cutree(fit, k = 3)
posterior_weights_df$id <- as.factor(1:nrow(posterior_weights_df))
posterior_weights_df$cluster <- clusters
melted_data <- melt(posterior_weights_df, id.vars = c("id", "cluster"))
melted_data$variable2 <- sub("_.*", "", melted_data$variable)
melted_data$variable3 <- as.factor(round(as.numeric(sub("*._", "", melted_data$variable)), 3))
melted_data <- melted_data %>% arrange(cluster)
ggplot(melted_data, aes(x = id, y = value, fill = variable3)) +
  geom_bar(stat = "identity") +
  facet_wrap(~cluster, scales = "free_y") +  # Facet by cluster
  labs(x = "ID", y = "Weight", fill = "PSD") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for better readability
  ggtitle("Structure Plot of Posterior Weights by Cluster") +
  coord_flip()  # Optional based on your preference for layout
ggsave(paste0(fig_dir, "/posterior_weights_structure_plot_clustered.pdf"), height = 10, width = 10)

### One miss classification: observation 18.
plot(datasets[[39]]$x, datasets[[39]]$y, type = 'p',  # Change to 'p' for points
     main = paste("Dataset", 39), xlab = "x", ylab = "y",
     cex = 0.3,
     cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
lines(datasets[[39]]$x, datasets[[39]]$truef, col = "red", lty = "dashed", lwd = 1.5)
lines(fitted_datasets[[39]]$x, fitted_datasets[[39]]$mean, col = "blue", lwd = 1.2)

plot(datasets[[18]]$x, datasets[[18]]$y, type = 'p',  # Change to 'p' for points
     main = paste("Dataset", 18), xlab = "x", ylab = "y",
     cex = 0.3,
     cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
lines(datasets[[18]]$x, datasets[[18]]$truef, col = "red", lty = "dashed", lwd = 1.5)
lines(fitted_datasets[[18]]$x, fitted_datasets[[18]]$mean, col = "blue", lwd = 1.2)

# visualize them:
pdf(paste0(fig_dir, "/fitted_datasets_catA.pdf"), height = 10, width = 10)
par(mfrow = c(4, 4), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (i in 1:16) {
  plot(datasets[[i]]$x, datasets[[i]]$y, type = 'p',  # Change to 'p' for points
       main = paste("Dataset", i), xlab = "x", ylab = "y",
       cex = 0.1,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  lines(datasets[[i]]$x, datasets[[i]]$truef, col = "red", lty = "dashed", lwd = 1.5)
  lines(fitted_datasets[[i]]$x, fitted_datasets[[i]]$mean, col = "blue", lwd = 1.2)
}
dev.off()
pdf(paste0(fig_dir, "/fitted_datasets_catB.pdf"), height = 10, width = 10)
par(mfrow = c(4, 4), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (i in 17:32) {
  plot(datasets[[i]]$x, datasets[[i]]$y, type = 'p',  # Change to 'p' for points
       main = paste("Dataset", i), xlab = "x", ylab = "y",
       cex = 0.1,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  lines(datasets[[i]]$x, datasets[[i]]$truef, col = "red", lty = "dashed", lwd = 1.5)
  lines(fitted_datasets[[i]]$x, fitted_datasets[[i]]$mean, col = "blue", lwd = 1.2)
}
dev.off()
pdf(paste0(fig_dir, "/fitted_datasets_catC.pdf"), height = 10, width = 10)
par(mfrow = c(4, 4), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (i in 33:48) {
  plot(datasets[[i]]$x, datasets[[i]]$y, type = 'p',  # Change to 'p' for points
       main = paste("Dataset", i), xlab = "x", ylab = "y",
       cex = 0.1,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  lines(datasets[[i]]$x, datasets[[i]]$truef, col = "red", lty = "dashed", lwd = 1.5)
  lines(fitted_datasets[[i]]$x, fitted_datasets[[i]]$mean, col = "blue", lwd = 1.2)
}
dev.off()
par(mfrow = c(1, 1))





