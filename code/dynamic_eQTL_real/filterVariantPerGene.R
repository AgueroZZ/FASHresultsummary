result_dir <- paste0(getwd(), "/output/dynamic_eQTL_real")
load(paste0(result_dir, "/fash_fit1_all.RData"))


datasets <- fash_fit1$fash_data$data_list
for (i in 1:length(datasets)) {
  datasets[[i]]$SE <- fash_fit1$fash_data$S[[i]]
}
all_genes <- unique(sapply(strsplit(names(datasets), "_"), "[[", 1))

all_genes
all_genes_vector <- sapply(strsplit(names(datasets), "_"), "[[", 1)
datasets_selected <- list()

compute_max_z_score <- function(index){
  dataset <- datasets[[index]]
  y <- dataset$y
  SE <- dataset$SE
  max(abs(y/SE))
}

for (gene in all_genes) {
  indices <- which(all_genes_vector == gene)
  results <- unlist(lapply(indices, compute_max_z_score))
  selected_index <- indices[which.max(results)]
  datasets_selected[length(datasets_selected) + 1] <- datasets[selected_index]
  names(datasets_selected)[length(datasets_selected)] <- names(datasets)[selected_index]
}

save(datasets_selected, file = "dynamicQTLs.RData")





