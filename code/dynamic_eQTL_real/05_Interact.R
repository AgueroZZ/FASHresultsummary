library(VariantAnnotation)
library(parallel)
num_cores <- 1

# find variant-gene pairs
dynamic_eqtl_sumstats <- read.delim("./data/dynamic_eqtl_sumstats.txt")
PC_data <- read.delim("./data/cell_line_ignore_missing_principal_components_9.txt")
rownames(PC_data) <- PC_data$Sample_id
quantile_normalized_no_proj <- read.csv("./data/quantile_normalized_no_projection.txt", sep="")
vcf <- readVcf("./data/YRI_genotype.vcf.gz", genome = "hg37")
geno_data <- geno(vcf)$DS
all_pairs <- dynamic_eqtl_sumstats[, c("variant_id", "gene_id")]
all_unique_genes <- unique(all_pairs$gene_id)

# A function that takes a gene_id, and analyze all its variants
analyze_one_gene_int <- function(selected_gene){
  all_variants_list <- list()
  all_relevant_variants <- all_pairs[all_pairs$gene_id == selected_gene, "variant_id"]
  all_time <- 0:15
  Y <- quantile_normalized_no_proj[which(quantile_normalized_no_proj$Gene_id == selected_gene),]

  linear_result <- data.frame()
  nonlinear_result <- data.frame()

  for (selected_variant in all_relevant_variants) {
    data_variant <- data.frame()
    for (selected_time in all_time) {
      selected_colnames <- grep(paste0("_", selected_time, "$"), colnames(Y), value = TRUE)
      Yt <- as.numeric(Y[,selected_colnames])
      selected_cell_lines <- gsub(paste0("^X|_", selected_time), "", selected_colnames)
      PC_selected <- PC_data[paste0(selected_cell_lines), c("PC1", "PC2", "PC3", "PC4", "PC5")]
      G <- geno_data[selected_variant, selected_cell_lines]
      # construct a data frame
      data_new <- data.frame(Yt, G, PC_selected)
      data_new$time <- selected_time
      data_variant <- rbind(data_variant, data_new)
    }
    data_variant$Yt <- (data_variant$Yt - mean(data_variant$Yt)) / sd(data_variant$Yt)
    # Test the non-linear results
    mod_fit <- lm(Yt ~ G*time + PC1*time + PC2*time + PC3*time + PC4*time + PC5*time +
                    G*I(time^2) + PC1*I(time^2) + PC2*I(time^2) + PC3*I(time^2) + PC4*I(time^2) + PC5*I(time^2)
                    , data = data_variant)

    mod_reduced <- lm(Yt ~ G + PC1*time + PC2*time + PC3*time + PC4*time + PC5*time +
                        PC1*I(time^2) + PC2*I(time^2) + PC3*I(time^2) + PC4*I(time^2) + PC5*I(time^2)
                      , data = data_variant)
    test_result <- lmtest::lrtest(mod_reduced, mod_fit)
    nonlinear_result_new <- data.frame(gene = selected_gene, variant = selected_variant,
                         pval = test_result$`Pr(>Chisq)`[2],
                         chisq = test_result$Chisq[2])

    data_variant$permuted_time <- sample(data_variant$time)
    mod_fit_permut <- lm(Yt ~ G + G:permuted_time + PC1*time + PC2*time + PC3*time + PC4*time + PC5*time +
                           G:I(permuted_time^2) + PC1*I(time^2) + PC2*I(time^2) + PC3*I(time^2) + PC4*I(time^2) + PC5*I(time^2)
                         , data = data_variant)
    test_result_perm <- lmtest::lrtest(mod_reduced, mod_fit_permut)
    nonlinear_result_new$permuted_pval <- test_result_perm$`Pr(>Chisq)`[2]
    nonlinear_result_new$permuted_chisq <- test_result_perm$Chisq[2]
    nonlinear_result <- rbind(nonlinear_result, nonlinear_result_new)

    # Test the linear results
    mod_fit <- lm(Yt ~ G*time + PC1*time + PC2*time + PC3*time + PC4*time + PC5*time, data = data_variant)
    mod_reduced <- lm(Yt ~ G + PC1*time + PC2*time + PC3*time + PC4*time + PC5*time, data = data_variant)
    test_result <- lmtest::lrtest(mod_reduced, mod_fit)
    linear_result_new <- data.frame(gene = selected_gene, variant = selected_variant,
                                       pval = test_result$`Pr(>Chisq)`[2],
                                       chisq = test_result$Chisq[2])
    mod_fit_permut <- lm(Yt ~ G + G:permuted_time + PC1*time + PC2*time + PC3*time + PC4*time + PC5*time, data = data_variant)
    test_result_perm <- lmtest::lrtest(mod_reduced, mod_fit_permut)
    linear_result_new$permuted_pval <- test_result_perm$`Pr(>Chisq)`[2]
    linear_result_new$permuted_chisq <- test_result_perm$Chisq[2]
    linear_result <- rbind(linear_result, linear_result_new)
  }
  return(list(linear_result = linear_result, nonlinear_result = nonlinear_result))
}

# Check the number of core to use
if(num_cores == 1){
  # Create a progress bar
  pb <- txtProgressBar(min = 0, max = length(all_unique_genes), style = 3)
  for (selected_gene in all_unique_genes) {
    setTxtProgressBar(pb, which(all_unique_genes == selected_gene))
    cat("Analyzing interactions in gene ", selected_gene, "\n")
    # check if the result already exists
    if (file.exists(paste0("./int_results/", selected_gene, ".rds"))) {
      cat("Result already exists, skipping gene ", selected_gene, "\n")
      next
    }
    result <- analyze_one_gene_int(selected_gene)
    saveRDS(result, paste0("./int_results/", selected_gene, ".rds"))
  }
}else{
  # use parallel package
  mclapply(all_unique_genes, function(selected_gene) {
    if (!file.exists(paste0("./int_results/", selected_gene, ".rds"))) {
      cat("Analyzing interactions in gene ", selected_gene, "\n")
      result <- analyze_one_gene_int(selected_gene)
      saveRDS(result, paste0("./int_results/", selected_gene, ".rds"))
    }else{
      cat("Result already exists, skipping gene ", selected_gene, "\n")
    }
  }, mc.cores = num_cores)
}



