library(tidyverse)

simulate_BF <- function(n = 1, from = 0, fit = 1, sebetahat = 1){
  # simulate a normal
  x <- rnorm(n = n, mean = 0, sd = from)

  # simulate the noise
  y <- x + rnorm(n, mean = 0, sd = sebetahat)

  # generate likelihood
  data <- ashr::set_data(betahat = y, sebetahat = sebetahat)
  ash_fit <- ashr::ash(betahat = y, sebetahat = sebetahat, mixcompdist = "normal", mixsd = fit, nullweight = 1, outputlevel = 3)
  ash_fit$fit_details$matrix_lik[,2]/ash_fit$fit_details$matrix_lik[,1]
}

BF0 <- simulate_BF(n = 8000, from = 0, fit = 0.8)
BF1 <- simulate_BF(n = 2000, from = 0.8, fit = 0.8)

## E_0(BF) = E_1(BF^{-1}) = 1
mean(BF0);mean(BF1^(-1))
## E_0(BF^k) = E_1(BF^{k-1})
mean(BF0^(-1)); mean(BF1^(-2))
mean(BF0^(2)); mean(BF1^(1))



## My hypothesis BF^2 does not have finite moment (indeed, only if fit > sebetahat)

n_vec <- c(100,500,1000,3000,10000)
B <- 100
sample_all <- NULL
for (n in n_vec){
  samples_n <- c()
  for (b in 1:B){
    BF <- simulate_BF(n = n, from = 0, fit = 0.8)
    mean_BFsquared <- mean(BF^(2))
    samples_n <- c(samples_n, mean_BFsquared)
  }
  sample_all <- cbind(sample_all, samples_n)
}
colnames(sample_all) <- n_vec

sample_all %>% apply(2, mean)
sample_all %>% apply(2, var)

1/sqrt(1 - (0.8^4))

# draw a boxplot for each column and put them side by side
boxplot((sample_all), col = "blue", xlab = "n", outline = F)

BF <- c(BF0, BF1)

(mean(BF) - 1.301448)/(1-1.301448)


# plot the histogram of BF0, BF1 and BF
par(mfrow = c(3,1))
hist(BF0, breaks = 100, col = "blue", xlim = c(0,10), xlab = "BF", main = "Histogram of BF0", freq = F)
hist(BF1, breaks = 100, col = "red", xlim = c(0,10), xlab = "BF", main = "Histogram of BF1", freq = F)
hist(BF, breaks = 100, col = "green", xlim = c(0,10), xlab = "BF", main = "Histogram of BF", freq = F)
par(mfrow = c(1,1))

par(mfrow = c(3,1))
hist(BF0^(-1), breaks = 100, col = "blue", xlim = c(0,2), xlab = "BF", main = "Histogram of BF0^(-1)", freq = F)
hist(BF1^(-1), breaks = 100, col = "red", xlim = c(0,2), xlab = "BF", main = "Histogram of BF1^(-1)", freq = F)
hist(BF^(-1), breaks = 100, col = "green", xlim = c(0,2), xlab = "BF", main = "Histogram of BF^(-1)", freq = F)
par(mfrow = c(1,1))

# plot the density of BF0, BF1 and BF
par(mfrow = c(3,1))
plot(density(BF0), xlim = c(0,10), col = "blue", xlab = "BF", main = "Density of BF0")
plot(density(BF1), xlim = c(0,10), col = "red", xlab = "BF", main = "Density of BF1")
plot(density(BF), xlim = c(0,10), col = "green", xlab = "BF", main = "Density of BF")
par(mfrow = c(1,1))



# A function to compute the running mean and record the index when it first touches 1

running_mean <- function(x, plot = T, cutoff = NULL){
  x <- sort(x, decreasing = F)

  if(is.null(cutoff)){
  n <- length(x)
  mean_x <- numeric(n)
  index <- numeric(n)
  for (i in 1:n){
    mean_x[i] <- mean(x[1:i])
    if (mean_x[i] >= 1){
      index[i] <- i
    }
  }
  if (plot){
    plot(mean_x, type = "l", xlab = "Index", ylab = "Running Mean", main = "Running Mean of BF")
    abline(h = 1, col = "red")
    abline(v = index[which(index != 0)][1], col = "blue")
  }
  return(list(mean = mean_x, index = index[which(index != 0)][1], cutoff = x[index[which(index != 0)][1]]))
  }
  else{
    return(list(mean = mean(x[x <= cutoff]), index = max(which(x <= cutoff)), cutoff = cutoff))
  }
}

cutoff_BF0 <- running_mean(BF)$cutoff
cutoff_BF1 <- running_mean(BF^(-1))$cutoff

running_mean(BF^(-1), cutoff = 1/cutoff_BF0)
running_mean(BF^(-1), cutoff = cutoff_BF1)

cutoff_H0 <- seq(0.1,10,0.1)
curve_1 <- c()
curve_2 <- c()
for (cutoff in cutoff_H0) {
  curve_1 <- c(curve_1, running_mean(BF, cutoff = cutoff)$mean)
  curve_2 <- c(curve_2, running_mean(BF^(-1), cutoff = cutoff^(-1))$mean)
}
plot(cutoff_H0, curve_1, type = "l", xlab = "Cutoff", ylab = "Running Mean", main = "Running Mean of BF", ylim = c(0,2))
lines(cutoff_H0, curve_2, col = "red")


# Region 1:
index_region1 <- which(BF <= (1/cutoff_BF1))
mean(index_region1 <= length(BF0))
mean(1:length(BF0) %in% index_region1)
# Region 2:
index_region2 <- which(BF > (1/cutoff_BF1) & BF <= cutoff_BF0)
mean(index_region2 <= length(BF0))
mean(1:length(BF0) %in% index_region2)
# Region 3:
index_region3 <- which(BF > cutoff_BF0)
mean(index_region3 <= length(BF0))
mean(1:length(BF0) %in% index_region3)



length(index_region1)/(length(index_region1) + length(index_region3))
# for each observation in region 2, compute a weight proportional to the distance to the boundary
distance_to_boundary1 <- cutoff_BF0 - BF[index_region2]
distance_to_boundary2 <- BF[index_region2] - (1/cutoff_BF1)
weight <- distance_to_boundary1/(distance_to_boundary1 + distance_to_boundary2)

(length(index_region1) + sum(weight))/(length(index_region1) + length(index_region2) + length(index_region3))
(length(index_region2) + 1 - sum(weight))/(length(index_region1) + length(index_region2) + length(index_region3))


# Compute the running mean for BF and BF^-1 based on the weight
sum(c(BF[index_region1], BF[index_region2]*weight))/(length(index_region1) + sum(weight))
sum(c((BF[index_region3]^(-1)), (BF[index_region2]^(-1))*(1-weight)))/(length(index_region3) + length(index_region2) - sum(weight))



compute_weighted_result <- function(BF, weight_scale = 1){
  cutoff_BF0 <- running_mean(BF, plot = FALSE)$cutoff
  cutoff_BF1 <- running_mean(BF^(-1), plot = FALSE)$cutoff

  index_region1 <- which(BF <= (1/cutoff_BF1))
  index_region2 <- which(BF > (1/cutoff_BF1) & BF <= cutoff_BF0)
  index_region3 <- which(BF > cutoff_BF0)

  distance_to_boundary1 <- weight_scale*(cutoff_BF0 - BF[index_region2])
  distance_to_boundary2 <- BF[index_region2] - (1/cutoff_BF1)
  weight <- distance_to_boundary1/(distance_to_boundary1 + distance_to_boundary2)

  return(
    list(
    BF0 = sum(c(BF[index_region1], BF[index_region2] * weight)) / (length(index_region1) + sum(weight)),
    invBF1 =  sum(c((BF[index_region3] ^ (-1)), (BF[index_region2] ^ (-1)) * (1 - weight))) / (length(index_region3) + length(index_region2) - sum(weight))
  )
  )
}

weight_vec <- exp(seq(-10,3,by = 1))
resultBF0 <- numeric(length(weight_vec))
resultBF1 <- numeric(length(weight_vec))
for (weight in weight_vec) {
  result <- compute_weighted_result(BF, weight_scale = weight)
  resultBF0[which(weight_vec == weight)] <- result$BF0
  resultBF1[which(weight_vec == weight)] <- result$invBF1
}

plot(weight_vec, resultBF0, type = "l",
     ylim = c(0,2),
     xlab = "Weight", ylab = "Weighted Mean", main = "Weighted Mean of BF0")
lines(weight_vec, resultBF1, col = "red")

diff <- 2 - resultBF0 - resultBF1

plot(weight_vec, diff, type = "l",
     xlab = "Weight", ylab = "Difference", main = "Difference between BF0 and BF1",
     ylim = c(-1,1))

running_mean(BF, cutoff = 1/cutoff_BF1)
running_mean(BF, cutoff = cutoff_BF0)

