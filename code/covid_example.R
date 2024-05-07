library(BayesGP)
library(TMB)
library(Matrix)
library(splines)
library(parallel)
library(reshape2)
library(mixsqp)
library(ISOweek)
library(lubridate)
library(tidyverse)

data_dir <- paste0(dirname(getwd()), "/data")
cpp_dir <- paste0(getwd(), "/cpp")
figure_dir <- paste0(dirname(getwd()), "/output/example")
result_dir <- paste0(dirname(getwd()), "/output/example/figure")
function_dir <- paste0(getwd(), "/function")

source(paste0(function_dir, "/functions_fitting_Poisson_covid.R"))
compile(paste0(cpp_dir, "/Poisson_covid.cpp"))
compile(paste0(cpp_dir, "/Poisson_just_fixed_covid.cpp"))
dyn.load(TMB::dynlib(paste0(cpp_dir, "/Poisson_covid")))
dyn.load(TMB::dynlib(paste0(cpp_dir, "/Poisson_just_fixed_covid")))

full_data_covid <- read.csv(file = paste0(data_dir, "/owid-covid-data.csv"), header = T)
full_data_covid <- full_data_covid %>% filter(date > "2020-01-01")
full_data_covid$y <- round(full_data_covid$new_deaths_per_million)
full_data_covid$quality <- abs((full_data_covid$y-full_data_covid$new_deaths_smoothed_per_million)/(full_data_covid$new_deaths_smoothed_per_million + 1))
# full_data_covid$y <- round(full_data_covid$new_deaths)
full_data_covid$Date <- as.Date(full_data_covid$date)
full_data_covid$x <- (as.Date(full_data_covid$date) %>% as.numeric())/31
## Assume COVID death rate is approximately 0 at time 2020-01-01, so set intercept being -3 
## Also, assume the derivatives are all zero at this point.
full_data_covid$x <- full_data_covid$x - (as.numeric(as.Date("2020-01-01"))/31)
# full_data_covid$x <- full_data_covid$x - min(full_data_covid$x)
full_data_covid$weekdays <- weekdays(as.Date(full_data_covid$date))
full_data_covid$weekdays <- factor(full_data_covid$weekdays,
                                    levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"),
                                    ordered = F)

## Create some extra_data_EU for date 2019-12-01 with zero death count:
full_data_covid_EU <- full_data_covid %>% filter(continent == "Europe") %>% select(Date, x, y, weekdays, location, iso_code, quality)
# extra_data_EU <- data.frame(Date = as.Date("2019-12-01"), x = 0, y = 0, weekdays = "Sunday", location = unique(full_data_covid_EU$location), iso_code = unique(full_data_covid_EU$iso_code), quality = 0)
# full_data_covid_EU <- rbind(full_data_covid_EU, extra_data_EU)

# Number of cores to use
num_cores <- 8

# Load the data
datasets <- list()
available_countries <- unique(full_data_covid_EU$location)
selected_countries <- c()
for (country in available_countries) {
  full_data <- full_data_covid_EU %>% filter(location == country, !is.na(y), !is.na(quality))
  full_data <- full_data %>% filter(quality <= 1.0, y != 0)
  if(sum(full_data$y != 0) < 300) {
    next
  }
  else{
    selected_countries <- c(selected_countries, country)
    datasets[[country]] <- full_data %>% arrange(Date)
  }
}

# visualize them:
par(mfrow = c(6, 6), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (i in 1:34) {
  plot(datasets[[i]]$x, (datasets[[i]]$y), type = 'p',  # Change to 'p' for points
       main = paste0(selected_countries[i]), xlab = "x", ylab = "y",
       cex = 0.1,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


# Compute L_vecs for each dataset in parallel
p_vec <- 3
log_prec_vec <- sort(unique(c(seq(0, 5, by = 0.1), seq(5, 10, by = 1), seq(-5,0, by = 1))))
# psd_iwp_vec <- seq(0.1,2, by = 0.1)
psd_iwp_vec <- exp(-.5*log_prec_vec)
L_vecs <- mclapply(datasets, function(dataset) {
  compute_log_likelihood_ospline_seq2(
    dataset = dataset,
    p = p_vec,
    num_knots = 50,
    psd_iwp_vector = psd_iwp_vec,
    pred_step = 1,
    betaprec = 0.0001
  )
}, mc.cores = num_cores)
L_matrix <- do.call(rbind, L_vecs)
max_indices <- apply(L_matrix, 1, which.max)
save(L_matrix, file = paste0(result_dir, "/L_matrix.rda"))
load(paste0(result_dir, "/L_matrix.rda"))

### Let's learn the weights through EB:
fit.sqp <- mixsqp(L = L_matrix, log = TRUE)
numiter <- nrow(fit.sqp$progress)
plot(1:numiter,fit.sqp$progress$objective,type = "b",
     pch = 20,lwd = 2,xlab = "SQP iteration",
     ylab = "objective",xaxp = c(1,numiter,numiter - 1))
fit.sqp$x ## the optimal weights
prior_weight <- data.frame(p = rep(p_vec, each = length(psd_iwp_vec)), psd_iwp = psd_iwp_vec, prior_weight = fit.sqp$x)

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
    dataset = dataset,
    num_knots = 50,
    prior_weight = prior_weight,
    betaprec = 0.0001,
    pred_step = 1
  )
  # Assume posterior weights are in the same order as in 'prior_weight$psd_iwp'
  posterior_weights_matrix[i, ] <- fit_result_final$posterior_weights[, "posterior_weight"]
  fitted_datasets[[i]] <- aggregate_fit_with_prior(x = dataset$x, fit_results_with_prior = fit_result_final, original = TRUE)$summary_df
}
names(fitted_datasets) <- selected_countries
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
melted_data$variable3 <- sub("*._", "", melted_data$variable)
melted_data$id <- selected_countries

ggplot(melted_data, aes(x = as.factor(id), y = value, fill = variable3)) +
  geom_bar(stat = "identity") +
  labs(x = "Country", y = "Weight", fill = "PSD") +
  theme_minimal() +
  ggtitle("Structure Plot of Posterior Weights") +
  coord_flip()  
ggsave(filename = paste0(figure_dir, "/structure_unclustered.pdf"), width = 10, height = 10)

# Normalize the data (optional but recommended)
wide_data <- scale(posterior_weights_matrix)
# Perform hierarchical clustering
d <- dist(wide_data, method = "euclidean")  # Distance matrix
fit <- hclust(d, method = "ward.D2")  # Clustering
clusters <- cutree(fit, k = 4)
# clusters <- as.numeric(factor(clusters, levels = c(2, 5, 3, 1, 6, 4), labels = c(1, 2, 3, 4, 5, 6)))

## Recode the factor so cluster 2 be 1, 1 be 2, 4 be 3, 3 be 4
clusters <- as.numeric(factor(clusters, levels = c(2, 1, 4, 3), labels = c(1, 2, 3, 4)))
melted_data$cluster <- clusters
posterior_weights_df$id <- 1:nrow(posterior_weights_df)
posterior_weights_df$cluster <- clusters
melted_data <- melt(posterior_weights_df, id.vars = c("id", "cluster"))
melted_data$variable2 <- sub("_.*", "", melted_data$variable)
melted_data$variable3 <- as.factor(round(as.numeric(sub("*._", "", melted_data$variable)), 3))
melted_data$id <- selected_countries
melted_data <- melted_data %>% arrange(cluster)
ggplot(melted_data, aes(x = id, y = value, fill = variable3)) +
  geom_bar(stat = "identity") +
  facet_wrap(~cluster, scales = "free_y") +  # Facet by cluster
  labs(x = "Country", y = "Weight", fill = "PSD") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels for better readability
  ggtitle("Structure Plot of Posterior Weights by Cluster") +
  coord_flip()  # Optional based on your preference for layout
ggsave(filename = paste0(figure_dir, "/structure_plot_clustered.pdf"), height = 10, width = 10)


### Plot of mortality curve for cluster
cluster1_countries <- selected_countries[clusters == 1]
cluster2_countries <- selected_countries[clusters == 2]
cluster3_countries <- selected_countries[clusters == 3]
cluster4_countries <- selected_countries[clusters == 4]

pdf(file = paste0(figure_dir, "/raw_data_cluster1.pdf"), width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (country in cluster1_countries[1:4]) {
  plot(datasets[[country]]$x, (datasets[[country]]$y), type = 'p',  # Change to 'p' for points
       main = paste0(country), xlab = "x", ylab = "y",
       cex = 0.5,
       cex.main = 1.5, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()
pdf(file = paste0(figure_dir, "/fit_data_cluster1.pdf"), width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (country in cluster1_countries[1:4]) {
    # fit_result_final <- fit_ospline_with_prior2(num_cores = num_cores, dataset = datasets[[country]], num_knots = 50, prior_weight = prior_weight, betaprec = 0.0001, pred_step = 1)
  agg_result <- fitted_datasets[[country]]
  plot(datasets[[country]]$Date, agg_result$median, type = 'l',  # Change to 'p' for points
       main = paste0(country), xlab = "", ylab = "", col = "blue",
       cex = 0.5, ylim = c(0,max(datasets[[country]]$y)*1.2),
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  points(datasets[[country]]$Date, datasets[[country]]$y, col = "black", cex = 0.1)
  polygon(c(datasets[[country]]$Date, rev(datasets[[country]]$Date)),
          c(agg_result$lower, rev(agg_result$upper)),
          col = rgb(0.6, 0.8, 1, alpha = 0.3), border = NA)
}
dev.off()

pdf(file = paste0(figure_dir, "/raw_data_cluster2.pdf"), width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (country in cluster2_countries[1:4]) {
  plot(datasets[[country]]$x, (datasets[[country]]$y), type = 'p',  # Change to 'p' for points
       main = paste0(country), xlab = "x", ylab = "y",
       cex = 0.5,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()
pdf(file = paste0(figure_dir, "/fit_data_cluster2.pdf"), width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (country in cluster2_countries[1:4]) {
    # fit_result_final <- fit_ospline_with_prior2(num_cores = num_cores, dataset = datasets[[country]], num_knots = 50, prior_weight = prior_weight, betaprec = 0.0001, pred_step = 1)
  agg_result <- fitted_datasets[[country]]
  plot(datasets[[country]]$Date, agg_result$median, type = 'l',  # Change to 'p' for points
       main = paste0(country), xlab = "", ylab = "", col = "blue",
       cex = 0.5, ylim = c(0,max(datasets[[country]]$y)*1.2),
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  points(datasets[[country]]$Date, datasets[[country]]$y, col = "black", cex = 0.1)
  polygon(c(datasets[[country]]$Date, rev(datasets[[country]]$Date)),
          c(agg_result$lower, rev(agg_result$upper)),
          col = rgb(0.6, 0.8, 1, alpha = 0.3), border = NA)
}
dev.off()

pdf(file = paste0(figure_dir, "/raw_data_cluster3.pdf"), width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (country in cluster3_countries[1:4]) {
  plot(datasets[[country]]$x, (datasets[[country]]$y), type = 'p',  # Change to 'p' for points
       main = paste0(country), xlab = "x", ylab = "y",
       cex = 0.5,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()
pdf(file = paste0(figure_dir, "/fit_data_cluster3.pdf"), width = 10, height = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (country in cluster3_countries[1:4]) {
    # fit_result_final <- fit_ospline_with_prior2(num_cores = num_cores, dataset = datasets[[country]], num_knots = 50, prior_weight = prior_weight, betaprec = 0.0001, pred_step = 1)
  agg_result <- fitted_datasets[[country]]
  plot(datasets[[country]]$Date, agg_result$median, type = 'l',  # Change to 'p' for points
       main = paste0(country), xlab = "", ylab = "", col = "blue",
       cex = 0.5, ylim = c(0,max(datasets[[country]]$y)*1.2),
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  points(datasets[[country]]$Date, datasets[[country]]$y, col = "black", cex = 0.1)
  polygon(c(datasets[[country]]$Date, rev(datasets[[country]]$Date)),
          c(agg_result$lower, rev(agg_result$upper)),
          col = rgb(0.6, 0.8, 1, alpha = 0.3), border = NA)
}
dev.off()


pdf(file = paste0(figure_dir, "/raw_data_cluster4.pdf"), width = 10, height = 10)
par(mfrow = c(1, 1), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (country in cluster4_countries[1:1]) {
  plot(datasets[[country]]$x, (datasets[[country]]$y), type = 'p',  # Change to 'p' for points
       main = paste0(country), xlab = "x", ylab = "y",
       cex = 0.5,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
}
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
dev.off()
pdf(file = paste0(figure_dir, "/fit_data_cluster4.pdf"), width = 10, height = 10)
par(mfrow = c(1, 1), mar = c(2, 2, 1, 1))  # Adjust margins as needed
for (country in cluster4_countries[1:1]) {
    # fit_result_final <- fit_ospline_with_prior2(num_cores = num_cores, dataset = datasets[[country]], num_knots = 50, prior_weight = prior_weight, betaprec = 0.0001, pred_step = 1)
  agg_result <- fitted_datasets[[country]]
  plot(datasets[[country]]$Date, agg_result$median, type = 'l',  # Change to 'p' for points
       main = paste0(country), xlab = "", ylab = "", col = "blue",
       cex = 0.5, ylim = c(0,max(datasets[[country]]$y)*1.2),
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)  # Adjust text size as needed
  points(datasets[[country]]$Date, datasets[[country]]$y, col = "black", cex = 0.1)
  polygon(c(datasets[[country]]$Date, rev(datasets[[country]]$Date)),
          c(agg_result$lower, rev(agg_result$upper)),
          col = rgb(0.6, 0.8, 1, alpha = 0.3), border = NA)
}
dev.off()

















