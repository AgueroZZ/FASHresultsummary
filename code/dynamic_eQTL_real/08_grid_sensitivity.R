library(fashr)
num_cores <- 16
load("./results/datasets_corrected.RData")


### Original grid: log_prec <- seq(0,10, by = 0.2)

### Finer grid:
log_prec <- seq(0,10, by = 0.1)
finer_grid <- sort(c(0, exp(-0.5*log_prec)))

fash_fit2_finer <- fash(Y = "beta", smooth_var = "time", S = "SE", data_list = datasets,
                  num_basis = 20, order = 2, betaprec = 0,
                  pred_step = 1, penalty = 10, grid = finer_grid,
                  num_cores = num_cores, verbose = TRUE)
save(fash_fit2_finer, file = "./results/fash_fit2_finer_all.RData")

### Sparser grid:
log_prec <- seq(0,10, by = 0.5)
sparse_grid <- sort(c(0, exp(-0.5*log_prec)))

fash_fit2_sparse <- fash(Y = "beta", smooth_var = "time", S = "SE", data_list = datasets,
                  num_basis = 20, order = 2, betaprec = 0,
                  pred_step = 1, penalty = 10, grid = sparse_grid,
                  num_cores = num_cores, verbose = TRUE)
save(fash_fit2_sparse, file = "./results/fash_fit2_sparse_all.RData")
