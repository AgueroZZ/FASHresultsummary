library(fashr)
num_cores <- 16

# Load the result from Strober et al.
result_strober <- read.delim("./data/dynamic_eqtl_sumstats.txt")

# all_gene we could study
all_genes_files <- list.files(path = "./results/")
all_genes <- gsub(pattern = "^X|.rds", "", all_genes_files)

datasets <- list()
for (gene in all_genes) {
  dataset <- readRDS(paste0("./results/", gene, ".rds"))
  names(dataset) <- paste0(gene, "_", names(dataset))
  datasets <- append(datasets, dataset)
}

# cat the total number of genes and the total number of datasets
cat("Total number of genes: ", length(all_genes), "\n")
cat("Total number of datasets: ", length(datasets), "\n")

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



