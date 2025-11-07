### Extract some representative genes from the cluster 9 that belong to `HYPOXIA` and `UV_RESPONSE_DN`
### also `IL2_STAT5_SIGNALING` and `ANGIOGENESIS`
# HYPOXIA_genes_entrez <- m_t2g$entrez_gene[m_t2g$gs_name == "HALLMARK_HYPOXIA"]
# UV_RESPONSE_DN_genes_entrez <- m_t2g$entrez_gene[m_t2g$gs_name == "HALLMARK_UV_RESPONSE_DN"]
# IL2_STAT5_SIGNALING_genes_entrez <- m_t2g$entrez_gene[m_t2g$gs_name == "HALLMARK_IL2_STAT5_SIGNALING"]
ANGIOGENESIS_genes_entrez <- m_t2g$entrez_gene[m_t2g$gs_name == "HALLMARK_ANGIOGENESIS"]

cluster9_genes <- sorted_posterior_weights_df %>% filter(cluster == 9) %>% pull(id)

cluster9_genes_converted <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = cluster9_genes,
  mart = mart
)

cluster9_genes_entrez <- cluster9_genes_converted$entrezgene_id
# overlap_genes_entrez <- intersect(HYPOXIA_genes_entrez, cluster9_genes_entrez)
# overlap_genes_entrez <- intersect(UV_RESPONSE_DN_genes_entrez, cluster9_genes_entrez)
# overlap_genes_entrez <- intersect(IL2_STAT5_SIGNALING_genes_entrez, cluster9_genes_entrez)
overlap_genes_entrez <- intersect(ANGIOGENESIS_genes_entrez, cluster9_genes_entrez)
overlap_genes <- cluster9_genes_converted$ensembl_gene_id[cluster9_genes_converted$entrezgene_id %in% overlap_genes_entrez]

fitted_datasets <- list()
pb <- txtProgressBar(min = 0, max = length(overlap_genes), style = 3)
for (i in 1:length(overlap_genes)) {
  gene <- overlap_genes[i]
  setTxtProgressBar(pb, i)
  dataset <- datasets[[gene]]
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
names(fitted_datasets) <- overlap_genes

set.seed(1234)
sample_index <- sample(1:length(overlap_genes), 4, replace = FALSE)
par(mfrow = c(2, 2))
for (i in 1:4) {
  plot_gene_expression_log_just_f(gene = overlap_genes[sample_index[i]], log_size_vec = log_size_vec)
}
par(mfrow = c(1, 1))


### As a comparision, also extract some genes from the cluster 1 that belong to `BILE_ACID_METABOLISM`
BILE_ACID_METABOLISM_genes_entrez <- m_t2g$entrez_gene[m_t2g$gs_name == "HALLMARK_BILE_ACID_METABOLISM"]

cluster1_genes <- sorted_posterior_weights_df %>% filter(cluster == 1) %>% pull(id)

cluster1_genes_converted <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = cluster1_genes,
  mart = mart
)

cluster1_genes_entrez <- cluster1_genes_converted$entrezgene_id
overlap_genes_entrez <- intersect(BILE_ACID_METABOLISM_genes_entrez, cluster1_genes_entrez)
overlap_genes <- cluster1_genes_converted$ensembl_gene_id[cluster1_genes_converted$entrezgene_id %in% overlap_genes_entrez]

fitted_datasets <- list()
pb <- txtProgressBar(min = 0, max = length(overlap_genes), style = 3)
for (i in 1:length(overlap_genes)) {
  gene <- overlap_genes[i]
  setTxtProgressBar(pb, i)
  dataset <- datasets[[gene]]
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
names(fitted_datasets) <- overlap_genes

set.seed(1234)
sample_index <- sample(1:length(overlap_genes), 4, replace = FALSE)
par(mfrow = c(2, 2))
for (i in 1:4) {
  plot_gene_expression_log_just_f(gene = overlap_genes[sample_index[i]], log_size_vec = log_size_vec)
}
par(mfrow = c(1, 1))



### Or compare with cluster 3 and `MYC_TARGETS_V2`
MYC_TARGETS_V2_genes_entrez <- m_t2g$entrez_gene[m_t2g$gs_name == "HALLMARK_MYC_TARGETS_V2"]
cluster3_genes <- sorted_posterior_weights_df %>% filter(cluster == 3) %>% pull(id)
cluster3_genes_converted <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = cluster3_genes,
  mart = mart
)
cluster3_genes_entrez <- cluster3_genes_converted$entrezgene_id
overlap_genes_entrez <- intersect(MYC_TARGETS_V2_genes_entrez, cluster3_genes_entrez)
overlap_genes <- cluster3_genes_converted$ensembl_gene_id[cluster3_genes_converted$entrezgene_id %in% overlap_genes_entrez]

fitted_datasets <- list()
pb <- txtProgressBar(min = 0, max = length(overlap_genes), style = 3)
for (i in 1:length(overlap_genes)) {
  gene <- overlap_genes[i]
  setTxtProgressBar(pb, i)
  dataset <- datasets[[gene]]
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
names(fitted_datasets) <- overlap_genes

set.seed(1234)
sample_index <- sample(1:length(overlap_genes), 4, replace = FALSE)
par(mfrow = c(2, 2))
for (i in 1:4) {
  plot_gene_expression_log_just_f(gene = overlap_genes[sample_index[i]], log_size_vec = log_size_vec)
}
par(mfrow = c(1, 1))





### Or compare with cluster 4 and `E2F_TARGETS`
E2F_TARGETS_genes_entrez <- m_t2g$entrez_gene[m_t2g$gs_name == "HALLMARK_E2F_TARGETS"]
cluster4_genes <- sorted_posterior_weights_df %>% filter(cluster == 4) %>% pull(id)
cluster4_genes_converted <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = cluster4_genes,
  mart = mart
)
cluster4_genes_entrez <- cluster4_genes_converted$entrezgene_id
overlap_genes_entrez <- intersect(E2F_TARGETS_genes_entrez, cluster4_genes_entrez)
overlap_genes <- cluster4_genes_converted$ensembl_gene_id[cluster4_genes_converted$entrezgene_id %in% overlap_genes_entrez]

fitted_datasets <- list()
pb <- txtProgressBar(min = 0, max = length(overlap_genes), style = 3)
for (i in 1:length(overlap_genes)) {
  gene <- overlap_genes[i]
  setTxtProgressBar(pb, i)
  dataset <- datasets[[gene]]
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
names(fitted_datasets) <- overlap_genes

set.seed(123)
sample_index <- sample(1:length(overlap_genes), 4, replace = FALSE)
par(mfrow = c(2, 2))
for (i in 1:4) {
  plot_gene_expression_log_just_f(gene = overlap_genes[sample_index[i]], log_size_vec = log_size_vec)
}
par(mfrow = c(1, 1))


### Or compare with cluster 8 and `EPITHELIAL_MESENCHYMAL_TRANSITION`

EPITHELIAL_MESENCHYMAL_TRANSITION_genes_entrez <- m_t2g$entrez_gene[m_t2g$gs_name == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
cluster8_genes <- sorted_posterior_weights_df %>% filter(cluster == 8) %>% pull(id)
cluster8_genes_converted <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id"),
  values = cluster8_genes,
  mart = mart
)
cluster8_genes_entrez <- cluster8_genes_converted$entrezgene_id
overlap_genes_entrez <- intersect(EPITHELIAL_MESENCHYMAL_TRANSITION_genes_entrez, cluster8_genes_entrez)
overlap_genes <- cluster8_genes_converted$ensembl_gene_id[cluster8_genes_converted$entrezgene_id %in% overlap_genes_entrez]

fitted_datasets <- list()
pb <- txtProgressBar(min = 0, max = length(overlap_genes), style = 3)
for (i in 1:length(overlap_genes)) {
  gene <- overlap_genes[i]
  setTxtProgressBar(pb, i)
  dataset <- datasets[[gene]]
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
names(fitted_datasets) <- overlap_genes

set.seed(123)
sample_index <- sample(1:length(overlap_genes), 4, replace = FALSE)
par(mfrow = c(2, 2))
for (i in 1:4) {
  plot_gene_expression_log_just_f(gene = overlap_genes[sample_index[i]], log_size_vec = log_size_vec)
}
par(mfrow = c(1, 1))
