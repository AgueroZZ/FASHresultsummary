library(VariantAnnotation)
library(parallel)
num_cores <- 16

# find variant-gene pairs
dynamic_eqtl_sumstats <- read.delim("./data/dynamic_eqtl_sumstats.txt")
PC_data <- read.delim("./data/principal_components_10.txt")
rownames(PC_data) <- PC_data$Sample_id
quantile_normalized_no_proj <- read.csv("./data/quantile_normalized_no_projection.txt", sep="")
vcf <- readVcf("./data/YRI_genotype.vcf.gz", genome = "hg37")
geno_data <- geno(vcf)$DS
all_pairs <- dynamic_eqtl_sumstats[, c("variant_id", "gene_id")]
all_unique_genes <- unique(all_pairs$gene_id)

# A function that takes a gene_id, and analyze all its variants
analyze_one_gene <- function(selected_gene){
  all_variants_list <- list()
  all_relevant_variants <- all_pairs[all_pairs$gene_id == selected_gene, "variant_id"]
  all_time <- 0:15
  Y <- quantile_normalized_no_proj[which(quantile_normalized_no_proj$Gene_id == selected_gene),]
  for (selected_variant in all_relevant_variants) {
    variant_result <- data.frame(time = all_time, gene = selected_gene, variant = selected_variant)
    beta_vec <- c(); SE_vec <- c(); pvalue_vec <- c()
    for (selected_time in all_time) {
      selected_colnames <- grep(paste0("_", selected_time, "$"), colnames(Y), value = TRUE)
      Yt <- Y[,selected_colnames]
      selected_cell_lines <- gsub(paste0("^X|_", selected_time), "", selected_colnames)
      PC_selected <- PC_data[paste0(selected_cell_lines, "_", selected_time), c("PC1", "PC2", "PC3", "PC4", "PC5")]
      G <- geno_data[selected_variant, selected_cell_lines]
      mod_fit <- lm(as.numeric(Yt) ~ G + PC_selected$PC1 + PC_selected$PC2 + PC_selected$PC3)
      beta_vec <- c(beta_vec, coef(mod_fit)[2])
      SE_vec <- c(SE_vec, summary(mod_fit)$coefficients[2,2])
      pvalue_vec <- c(pvalue_vec, summary(mod_fit)$coefficients[2,4])
    }
    variant_result$beta <- beta_vec
    variant_result$SE <- SE_vec
    variant_result$pvalue <- pvalue_vec
    all_variants_list[[selected_variant]] <- variant_result
  }

  return(all_variants_list)
}

# Check the number of core to use
if(num_cores == 1){
  # Create a progress bar
  pb <- txtProgressBar(min = 0, max = length(all_unique_genes), style = 3)
  for (selected_gene in all_unique_genes) {
    setTxtProgressBar(pb, which(all_unique_genes == selected_gene))
    cat("Analyzing gene ", selected_gene, "\n")
    # check if the result already exists
    if (file.exists(paste0("./eQTL_results/", selected_gene, ".rds"))) {
      cat("Result already exists, skipping gene ", selected_gene, "\n")
      next
    }
    result <- analyze_one_gene(selected_gene)
    saveRDS(result, paste0("./eQTL_results/", selected_gene, ".rds"))
  }
}else{
  # use parallel package
  mclapply(all_unique_genes, function(selected_gene) {
    if (!file.exists(paste0("./eQTL_results/", selected_gene, ".rds"))) {
      cat("Analyzing gene ", selected_gene, "\n")
      result <- analyze_one_gene(selected_gene)
      saveRDS(result, paste0("./eQTL_results/", selected_gene, ".rds"))
    }else{
      cat("Result already exists, skipping gene ", selected_gene, "\n")
    }
  }, mc.cores = num_cores)
}


