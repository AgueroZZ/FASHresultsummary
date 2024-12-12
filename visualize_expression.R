library(BayesGP)
library(TMB)
library(Matrix)
library(splines)
library(parallel)
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

plot_gene_expression <- function(gene, cell_line_index) {
  load(paste0(result_dir, "/fitted_pos_results", cell_line_index, ".rda"))
  load(paste0(result_dir, "/datasets_", cell_line_index, ".rda"))
  agg_result <- fitted_datasets[[gene]]
  plot(datasets[[gene]]$Day, agg_result$median, type = 'l',
       main = paste0(gene), xlab = "Day", ylab = "Expression", col = "blue",
       cex = 0.5, ylim = c(0, max(datasets[[gene]]$y) * 1.2),
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
  points(datasets[[gene]]$Day, datasets[[gene]]$y, col = "black", cex = 0.5)
  polygon(c(datasets[[gene]]$Day, rev(datasets[[gene]]$Day)),
          c(agg_result$lower, rev(agg_result$upper)),
          col = rgb(0.6, 0.8, 1, alpha = 0.3), border = NA)
}
plot_gene_expression_across_cell_lines <- function(gene, cell_line_indices) {
  # Initialize an empty list to store data for plotting
  plot_data <- list()

  # Loop over each cell line index
  for (cell_line_index in cell_line_indices) {
    # Load results and datasets for the current cell line
    result_file <- paste0(result_dir, "/fitted_pos_results", cell_line_index, ".rda")
    dataset_file <- paste0(result_dir, "/datasets_", cell_line_index, ".rda")

    if (file.exists(result_file) && file.exists(dataset_file)) {
      load(result_file)
      load(dataset_file)

      # Extract aggregated results and days for the specified gene
      if (gene %in% names(fitted_datasets) && gene %in% names(datasets)) {
        agg_result <- fitted_datasets[[gene]]
        day <- datasets[[gene]]$Day
        plot_data[[cell_line_index]] <- data.frame(
          Day = day,
          Median = agg_result$median,
          CellLine = paste0("Cell Line ", cell_line_index)
        )
      } else {
        message(paste("Gene", gene, "not found in cell line", cell_line_index))
      }
    } else {
      message(paste("Data not found for cell line", cell_line_index))
    }
  }

  # Combine all data into a single data frame
  plot_data_combined <- bind_rows(plot_data)

  # Generate the plot using ggplot2
  ggplot(plot_data_combined, aes(x = Day, y = Median, color = CellLine, group = CellLine)) +
    geom_line(size = 1) +
    labs(
      title = paste("Posterior Median of", gene, "Across Cell Lines"),
      x = "Day",
      y = "Expression",
      color = "Cell Line"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
}
# plot_gene_expression_across_cell_lines("ENSG00000070371",cell_line_indices = c(1:19))

plot_gene_expression_log_just_f <- function(gene, cell_line_index) {
  load(paste0(result_dir, "/fitted_pos_results", cell_line_index, ".rda"))
  load(paste0(result_dir, "/datasets_", cell_line_index, ".rda"))
  size_vec <- numeric(length = nrow(datasets[[1]]))
  for (i in 1:length(size_vec)) {
    all_counts <- unlist(lapply(datasets, function(x) {x$y[x$Day == (i-1)]}))
    size_vec[i] <- sum(all_counts)
  }
  log_size_vec <- log(size_vec)
  agg_result <- fitted_datasets[[gene]]
  agg_result[,c(2,3,4)] <-  agg_result[,c(2,3,4)]/exp(log_size_vec)
  min_y <- min(agg_result$lower)
  max_y <- max(agg_result$upper)  # Assuming this is positive and makes sense in your context

  # Set the limits for the y-axis
  log_min_y <- min_y / 2
  log_max_y <- max_y * 1.2

  plot(datasets[[gene]]$Day, agg_result$median, type = 'l',
       main = paste0(gene), xlab = "Day", ylab = "Expression", col = "blue",
       cex = 0.5,
       cex.main = 0.8, cex.lab = 0.7, cex.axis = 0.7)
  polygon(c(datasets[[gene]]$Day, rev(datasets[[gene]]$Day)),
          c(agg_result$lower, rev(agg_result$upper)),
          col = rgb(0.6, 0.8, 1, alpha = 0.3), border = NA)
}
plot_gene_expression_log_across_cell_lines <- function(gene, cell_line_indices, color = TRUE, CI = TRUE) {
  # Initialize an empty list to store data for plotting
  plot_data <- list()

  # Loop over each cell line index
  for (cell_line_index in cell_line_indices) {
    # Load results and datasets for the current cell line
    result_file <- paste0(result_dir, "/fitted_pos_results", cell_line_index, ".rda")
    dataset_file <- paste0(result_dir, "/datasets_", cell_line_index, ".rda")

    if (file.exists(result_file) && file.exists(dataset_file)) {
      load(result_file)
      load(dataset_file)

      # Calculate size vector for normalization
      size_vec <- numeric(length = nrow(datasets[[1]]))
      for (i in seq_along(size_vec)) {
        all_counts <- unlist(lapply(datasets, function(x) { x$y[x$Day == (i - 1)] }))
        size_vec[i] <- sum(all_counts)
      }
      log_size_vec <- log(size_vec)

      # Extract aggregated results and normalize by size vector
      if (gene %in% names(fitted_datasets) && gene %in% names(datasets)) {
        agg_result <- fitted_datasets[[gene]]
        agg_result[, c("median", "lower", "upper")] <- agg_result[, c("median", "lower", "upper")] / exp(log_size_vec)

        # Prepare data for plotting
        plot_data[[cell_line_index]] <- data.frame(
          Day = datasets[[gene]]$Day,
          Median = agg_result$median,
          Lower = agg_result$lower,
          Upper = agg_result$upper,
          CellLine = paste0("Cell Line ", cell_line_index)
        )
      } else {
        message(paste("Gene", gene, "not found in cell line", cell_line_index))
      }
    } else {
      message(paste("Data not found for cell line", cell_line_index))
    }
  }

  # Combine all data into a single data frame
  plot_data_combined <- bind_rows(plot_data)

  # Base ggplot
  p <- ggplot(plot_data_combined, aes(x = Day, y = Median, group = CellLine)) +
    geom_line(size = 1, aes(color = if (color) CellLine else NULL)) +
    labs(
      title = paste("Normalized Log Expression of", gene, "Across Cell Lines"),
      x = "Day",
      y = "Normalized Log Expression"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )

  # Apply single pink color if color = FALSE
  if (!color) {
    pink_color <- rgb(1, 0.4, 0.6, alpha = 0.6)  # Define a pink color with transparency
    p <- p +
      geom_line(size = 1, color = pink_color) +  # Set line color explicitly
      guides(color = "none")
  }

  # Add confidence intervals if CI is TRUE
  if (CI) {
    if (color) {
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = CellLine), alpha = 0.2, color = NA)
    } else {
      pink_color <- rgb(1, 0.4, 0.6, alpha = 0.2)  # Transparent pink for ribbon
      p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = pink_color, alpha = 0.2, color = NA) +
        guides(fill = "none")
    }
  }

  # Return the plot
  return(p)
}
plot_gene_expression_log_across_cell_lines("ENSG00000070371", cell_line_indices = c(1:5), color = F, CI = F)

### Compute the posterior median of sigma for each gene
compute_sigma_means <- function(cell_line_indices) {

  # Take out the indices that are not in the range of cell_line_indices
  cell_line_indices <- cell_line_indices[cell_line_indices >= 1 & cell_line_indices <= 19]
  # Take out the indices whose results are not available
  cell_line_indices <- cell_line_indices[-c(7,8,10,14,18,19)]

  sigma_means_list <- list()  # Initialize list to store sigma_mean vectors

  # Loop over each cell line index
  for (cell_line_index in cell_line_indices) {
    # Construct the path to the result file
    result_file <- paste0(result_dir, "/L_matrix_", cell_line_index, ".rda")
    if (file.exists(result_file)) {
      # Load the result file
      load(result_file)
      print(cell_line_index)
      fit.sqp <- mixsqp(L = L_matrix, log = TRUE, control = list(verbose = FALSE))
      posterior_matrix <- matrix(0, nrow = nrow(L_matrix), ncol = ncol(L_matrix))
      for(i in 1:nrow(L_matrix)){
        posterior_matrix[i,] <- exp(L_matrix[i,] - max(L_matrix[i,]) + log(fit.sqp$x))
        posterior_matrix[i,] <- posterior_matrix[i,]/sum(posterior_matrix[i,])
      }
      colnames(posterior_matrix) <- paste0(2,"_",psd_vec)

      # Extract sigma_val from column names of posterior_matrix
      sigma_val <- colnames(posterior_matrix)
      sigma_val <- as.numeric(sub(".*_", "", sigma_val))

      # Compute sigma_mean for the current cell line
      sigma_mean <- posterior_matrix %*% sigma_val

      # Add sigma_mean to the list
      sigma_means_list[[cell_line_index]] <- sigma_mean
    } else {
      message(paste("Result file not found for cell line index:", cell_line_index))
    }
  }

  # Combine sigma_mean vectors into a matrix
  sigma_means_matrix <- do.call(cbind, sigma_means_list)
  colnames(sigma_means_matrix) <- paste0("CellLine_", cell_line_indices)
  rownames(sigma_means_matrix) <- all_gene

  return(sigma_means_matrix)
}
sigma_matrix <- compute_sigma_means(c(1:19))









