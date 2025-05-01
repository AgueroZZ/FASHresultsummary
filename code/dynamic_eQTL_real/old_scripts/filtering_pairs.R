# file path 
file_path <- "./data/strober_nondyn/"
# all files
all_files <- list.files(path = file_path)

# merge all files into one single data, by rbinding them
data <- do.call(rbind, lapply(all_files, function(x) read.delim(paste0(file_path, x))))


select_gene_variant <- function(data, pair_threshold = 1){
  selected_data <- data[data$pvalue <= pair_threshold, ]
  # select genes that we want to study
  selected_genes <- unique(selected_data$ensamble_id)
  # select variants for each gene
  selected_variants <- list()
  # create a progress bar
  pb <- txtProgressBar(min = 0, max = length(selected_genes), style = 3)
  for (gene in selected_genes) {
    # update progress bar
    setTxtProgressBar(pb, which(selected_genes == gene))
    new_variants <- list(unique(selected_data$rs_id[selected_data$ensamble_id == gene]))
    names(new_variants) <- gene
    selected_variants <- append(selected_variants, new_variants)
  }
  return(list(selected_genes = selected_genes, selected_variants = selected_variants))
}


selected_pairs <- select_gene_variant(data, pair_threshold = 0.05/16)

# sum the total number of variant and gene pairs
total_pairs_selected <- sum(sapply(selected_pairs$selected_variants, length))
total_pairs_selected/1009173

# total number of genes
length(selected_pairs$selected_genes)

# output the selected pairs
saveRDS(selected_pairs, file = "./results/selected_pairs.rds")