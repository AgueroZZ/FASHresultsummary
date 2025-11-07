library(fashr)

# Simulate five datasets
set.seed(123)
datasets <- list()
for (i in 1:200) {
  n <- 20
  t <- seq(0, 5, length.out = n)
  sd <- sample(c(0.5, 1), size = n, replace = TRUE)
  u <- runif(1); if (u < 0.6) {f <- function(t) 2*u*cos(t)} else {f <- function(t) (u)}
  y <- f(t) + rnorm(n, sd = sd)
  datasets[[i]] <- data.frame(t = t, y = y, sd = sd, category = ifelse(u >= 0.6, "H0", "H1"))
}

# Fit the model
fash_fit <- fash(Y = "y", smooth_var = "t", S = "sd", data = datasets,
                 betaprec = 0, penalty = 2, grid = c(0, 0.1),
                 order = 1, likelihood = "gaussian", verbose = TRUE)
fash_fit$prior_weights

FDR_result <- fdr_control(fash_fit)
answers <- unlist(lapply(datasets, function(x) x$category[1]))
discovers <- FDR_result$fdr_results$index[FDR_result$fdr_results$FDR < 0.1]
mean(answers[discovers] == "H0")

# try min-lfsr
t_vec <- seq(0, 5, length.out = 100)
min_lfsr <- numeric(length(datasets))
# create a progress bar
pb <- txtProgressBar(min = 0, max = length(datasets), style = 3)
for (i in 1:length(datasets)) {
  setTxtProgressBar(pb, i)
  selected_index <- i
  sample_i <- predict(fash_fit, index = selected_index, smooth_var = t_vec, only.samples = T)
  sample_i_tilde <- apply(sample_i, 2, function(x) x - x[1])
  pos_prob <- apply(sample_i_tilde, 1, function(x) mean(x >= 0))
  neg_prob <- apply(sample_i_tilde, 1, function(x) mean(x <= 0))
  lfsr <- pmin(pos_prob, neg_prob)
  min_lfsr[i] <- min(lfsr)
}
lfsr_df <- data.frame(min_lfsr = min_lfsr, index = 1:length(datasets))
lfsr_df <- lfsr_df[order(lfsr_df$min_lfsr),]
lfsr_df$fsr <- cumsum(lfsr_df$min_lfsr) / (1:length(datasets))
discovers_fsr <- lfsr_df$index[lfsr_df$fsr < 0.1]
mean(answers[discovers_fsr] == "H0")

# compare min-lfsr against lfdr
plot(min_lfsr, fash_fit$posterior_weights[,1],
     ylim = c(0, 1), xlim = c(0, 1),
     xlab = "min-lfsr", ylab = "lfdr")
abline(a = 0, b = 1, col = "red")


