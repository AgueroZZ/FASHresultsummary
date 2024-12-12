.libPaths(c("~/lib", .libPaths()))

library(BayesGP)
library(TMB)
library(Matrix)
library(splines)
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)
library(reshape2)
library(mixsqp)
library(tidyverse)

cpp_dir <- paste0(getwd(), "/code/cpp")
fig_dir <- paste0(getwd(), "/output/expression")
result_dir <- paste0(getwd(), "/output/expression")
function_dir <- paste0(getwd(), "/code/function")
data_dir <- paste0(getwd(), "/data/expression_data")
source(paste0(function_dir, "/functions_fitting_Poisson_expression.R"))
compile(paste0(cpp_dir, "/Poisson_expression.cpp"))
compile(paste0(cpp_dir, "/Poisson_just_fixed_expression.cpp"))
dyn.load(TMB::dynlib(paste0(cpp_dir, "/Poisson_expression")))
dyn.load(TMB::dynlib(paste0(cpp_dir, "/Poisson_just_fixed_expression")))
num_cores <- 4
betaprec <- 1e-6

all_data_names <- list.files(data_dir)

pipeline_once <- function(cell_line_index){
  # Step 1: Load the datasets
  load(paste0(data_dir, "/", all_data_names[cell_line_index]))
  num_knots <- 16
  p <- 2
  datasets <- list()
  all_gene <- unique(expression_data_one_cell_line$Gene_id)
  # Check if the output already exists
  if (!file.exists(paste0(result_dir, "/datasets_", cell_line_index, ".rda"))) {
    for (gene in all_gene) {
      datasets[[gene]] <- expression_data_one_cell_line %>% filter(Gene_id == gene)
      datasets[[gene]]$x <- as.numeric(datasets[[gene]]$Day)
      datasets[[gene]]$y <- datasets[[gene]]$value
    }
    save(datasets, file = paste0(result_dir, "/datasets_", cell_line_index, ".rda"))
  }else {
    load(paste0(result_dir, "/datasets_", cell_line_index, ".rda"))
  }

  # Step 2: Compute the library size
  size_vec <- numeric(length = nrow(datasets[[1]]))
  for (i in 1:length(size_vec)) {
    Day <- datasets[[1]]$Day[i]
    all_counts <- unlist(lapply(datasets, function(x) {x$y[x$Day == Day]}))
    size_vec[i] <- sum(all_counts)
  }
  log_size_vec <- log(size_vec)

  # Step 3: Performing the EB
  psd_vec <- sort(unique(c(0, seq(0,3, length.out = 31))))
  L_vecs <- vector("list", length(datasets))
  pb <- txtProgressBar(min = 0, max = length(datasets), style = 3)
  # Check if L_matrix already exists
  if (file.exists(paste0(result_dir, "/L_matrix_", cell_line_index, ".rda"))) {
    load(paste0(result_dir, "/L_matrix_", cell_line_index, ".rda"))
  }else{
    for (i in 1:length(datasets)) {
      setTxtProgressBar(pb, i)
      dataset <- datasets[[i]]
      L_vecs[[i]] <- compute_log_likelihood_ospline_seq2(
        dataset = dataset,
        p = p,
        num_knots = num_knots,
        psd_iwp_vector = psd_vec,
        pred_step = 1,
        betaprec = betaprec,
        log_lib_size = log_size_vec
      )
    }
    L_matrix <- do.call(rbind, L_vecs)
    save(L_matrix, file = paste0(result_dir, "/L_matrix_", cell_line_index, ".rda"))
  }

  # Step 4: Fitting the model
  fit.sqp <- mixsqp(L = L_matrix, log = TRUE, control = list(verbose = FALSE))
  numiter <- nrow(fit.sqp$progress)
  prior_weight <- data.frame(p = rep(p, each = length(psd_vec)), psd_iwp = psd_vec, prior_weight = fit.sqp$x)
  posterior_matrix <- matrix(0, nrow = nrow(L_matrix), ncol = ncol(L_matrix))
  for(i in 1:nrow(L_matrix)){
    posterior_matrix[i,] <- exp(L_matrix[i,] - max(L_matrix[i,]) + log(fit.sqp$x))
    posterior_matrix[i,] <- posterior_matrix[i,]/sum(posterior_matrix[i,])
  }
  colnames(posterior_matrix) <- paste0(p,"_",psd_vec)
  posterior_weights_matrix <- posterior_matrix
  num_datasets <- length(datasets)
  num_weights <- sum(prior_weight$prior_weight != 0)
  # Check if the output already exists
  if (!file.exists(paste0(result_dir, "/fitted_pos_results", cell_line_index, ".rda"))) {
    # Loop through each dataset and perform fitting
    fitted_datasets <- list()
    # create a progress bar
    pb <- txtProgressBar(min = 0, max = num_datasets, style = 3)
    for (i in seq_along(datasets)) {
      setTxtProgressBar(pb, i)
      dataset <- datasets[[i]]
      fit_result_final <- fit_ospline_with_prior2(
        num_cores = 1,
        dataset = dataset,
        num_knots = num_knots,
        prior_weight = prior_weight,
        betaprec = betaprec,
        pred_step = 1,
        log_lib_size = log_size_vec
      )
      fitted_datasets[[i]] <- aggregate_fit_with_prior(x = dataset$x, fit_results_with_prior = fit_result_final, original = TRUE)$summary_df
    }
    names(fitted_datasets) <- all_gene
    save(fitted_datasets, posterior_weights_matrix,
         file = paste0(result_dir, "/fitted_pos_results", cell_line_index, ".rda"))
  }else{
    load(paste0(result_dir, "/fitted_pos_results", cell_line_index, ".rda"))
  }
}


num_cores <- 8  # Reserve one core for system processes
cl <- makeCluster(num_cores)

makeCluster(num_cores)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths("~/lib"))

# Export required variables to each worker node
clusterExport(cl, c("pipeline_once", "cpp_dir", "result_dir", "data_dir",
                    "all_data_names", "betaprec"))

# Run the paralleled loop
foreach(cell_line_index = 1:length(all_data_names),
        .packages = c("BayesGP", "TMB", "Matrix", "splines", "parallel",
                      "doParallel", "foreach", "ggplot2", "reshape2",
                      "mixsqp", "tidyverse")) %dopar% {

                        dyn.load(TMB::dynlib(paste0(cpp_dir, "/Poisson_expression")))
                        dyn.load(TMB::dynlib(paste0(cpp_dir, "/Poisson_just_fixed_expression")))

                        pipeline_once(cell_line_index)
                      }

stopCluster(cl)
registerDoSEQ()
