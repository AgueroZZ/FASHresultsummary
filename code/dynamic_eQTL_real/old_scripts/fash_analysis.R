library(fashr)
num_cores <- 8

# Load the result from Strober et al.
result_strober <- read.delim("./data/dynamic_eqtl_sumstats.txt")

# all_gene we could study
all_genes_files <- list.files(path = "./results/")
all_genes <- gsub(pattern = "^X|.rds", "", all_genes_files)

# load the selected pairs
selected_pairs <- readRDS("./eQTL_results/selected_pairs.rds")
selected_genes <- selected_pairs$selected_genes

datasets <- list()
for (gene in selected_genes) {
  dataset <- readRDS(paste0("./results/", gene, ".rds"))
  
  # only keep the dataset in selected_pairs
  dataset_kept <- dataset[selected_pairs$selected_variants[[gene]]]
  
  # augment the dataset with the gene name
  names(dataset_kept) <- paste0(gene, "_", names(dataset_kept))
  datasets <- append(datasets, dataset_kept)
}

# cat the total number of genes and the total number of datasets
cat("Total number of genes: ", length(selected_genes), "\n")
cat("Total number of datasets: ", length(datasets), "\n")

log_prec <- seq(0,10, by = 0.2)
fine_grid <- sort(c(0, exp(-0.5*log_prec)))
# plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "Time", ylab = "Beta")
# # draw fine grid
# abline(v = fine_grid, col = "grey", lty = 2)


fash_fit1 <- fash(Y = "beta", smooth_var = "time", S = "SE", data_list = datasets,
                 num_basis = 20, order = 1, betaprec = 1e-6,
                 pred_step = 1, penalty = 30, grid = fine_grid,
                 num_cores = num_cores, verbose = TRUE)

save(fash_fit1, file = "./eQTL_results/fash_fit1.RData")

fash_fit2 <- fash(Y = "beta", smooth_var = "time", S = "SE", data_list = datasets,
                 num_basis = 20, order = 2, betaprec = 1e-6,
                 pred_step = 1, penalty = 30, grid = fine_grid,
                 num_cores = num_cores, verbose = TRUE)

save(fash_fit2, file = "./eQTL_results/fash_fit2.RData")






