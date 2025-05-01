library(fashr)
num_cores <- 16

# Load the result from Strober et al.
result_strober <- read.delim("./data/dynamic_eqtl_sumstats.txt")

# all_gene we could study
all_genes_files <- list.files(path = "./eQTL_results/")
all_genes <- gsub(pattern = "^X|.rds", "", all_genes_files)

datasets <- list()
for (gene in all_genes) {
  dataset <- readRDS(paste0("./eQTL_results/", gene, ".rds"))
  names(dataset) <- paste0(gene, "_", names(dataset))
  datasets <- append(datasets, dataset)
}

# cat the total number of genes and the total number of datasets
cat("Total number of genes: ", length(all_genes), "\n")
cat("Total number of datasets: ", length(datasets), "\n")


# Trying the t-based correction to SE:
num_vec <- c(19, 19, 16, 19, 16, 19, 19, 19, 19, 19, 19, 19, 19, 18, 19, 19)
# A function that takes a vector of y and a vector of SE, then correct the SE
correct_SE <- function(y_vec, SE_vec, df_vec){
  # p-value under t-distribution
  t_val <- y_vec / SE_vec
  p_val <- 2 * pt(abs(t_val), df = df_vec, lower.tail = FALSE)
  # z-score
  z_val <- qnorm(1 - p_val / 2)
  z_val <- z_val * sign(t_val)
  SE_vec_corrected <- abs(y_vec) / abs(z_val)
  SE_vec_corrected
}
for (i in 1:length(datasets)) {
  SE_vec_corrected <- correct_SE(y_vec = datasets[[i]]$beta,
                                 SE_vec = datasets[[i]]$SE,
                                 df_vec = (num_vec - 5))
  # update the dataset
  datasets[[i]]$SE <- SE_vec_corrected
}

save(datasets, file = "./results/datasets_corrected.RData")


load("./results/datasets_corrected.RData")

log_prec <- seq(0,10, by = 0.2)
fine_grid <- sort(c(0, exp(-0.5*log_prec)))
fine_grid

fash_fit1 <- fash(Y = "beta", smooth_var = "time", S = "SE", data_list = datasets,
                  num_basis = 20, order = 1, betaprec = 0,
                  pred_step = 1, penalty = 10, grid = fine_grid,
                  num_cores = num_cores, verbose = TRUE)
save(fash_fit1, file = "./results/fash_fit1_all.RData")

fash_fit2 <- fash(Y = "beta", smooth_var = "time", S = "SE", data_list = datasets,
                  num_basis = 20, order = 2, betaprec = 0,
                  pred_step = 1, penalty = 10, grid = fine_grid,
                  num_cores = num_cores, verbose = TRUE)
save(fash_fit2, file = "./results/fash_fit2_all.RData")



