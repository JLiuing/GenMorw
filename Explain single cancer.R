## Retrieval of contribution percentages from various PPIs to cancer cohort prioritization genes
# After completing the prioritization of a cancer cohort, the program can assist users in analyzing 
# the contribution percentage of different PPI networks to the gene prioritization. By specifying 
# the "cancer" parameter as the TCGA abbreviation for the cancer, this analysis can be performed.

cancer <- "LAML"

library(tidyverse)
source("AvailablePPI.R")
PPIs <- AvailablePPI()
results_of_PPIs <- list()
all_considered_genes <- c()
PPI_contribution <- list()
PPI_contribution_norm <- list()
for (PPI in PPIs)
{
  patient_results_symbol <- readRDS(paste0("Results/Data/", cancer, "/", PPI, "/patient_results_symbol.rds"))
  all_considered_genes_PPI <- c()
  gene_sorts <- list()
  for (patient in patient_results_symbol)
  {
    if (length(patient) == 0) {
      next
    }
    tsort <- 1:length(patient)
    tsort <- (tsort / length(patient)) * log2(1 + c(1:length(patient)))
    dataFrame <- data.frame(gene = patient, sort = tsort)
    rownames(dataFrame) <- dataFrame$gene
    for (mutgene in patient)
    {
      gene_sorts[[mutgene]] <- c(gene_sorts[[mutgene]], dataFrame[mutgene, 2])
    }
    all_considered_genes_PPI <- c(all_considered_genes_PPI, patient)
  }
  rm(patient, dataFrame, mutgene)

  all_considered_genes_PPI <- all_considered_genes_PPI %>%
    unique() %>%
    sort()
  final_sort <- c()
  for (gene in all_considered_genes_PPI)
  {
    score <- mean(gene_sorts[[gene]]) / length(gene_sorts[[gene]])
    final_sort <- c(final_sort, score)
  }
  rm(gene, score)
  names(final_sort) <- all_considered_genes_PPI
  final_sort <- final_sort %>% sort()
  results_of_PPIs[[PPI]] <- final_sort %>% names()
  all_considered_genes <- c(all_considered_genes, results_of_PPIs[[PPI]])

  tp_PPI_contribution <- list()
  tp_PPI_contribution_norm <- list()
  for (gene in all_considered_genes_PPI)
  {
    tp <- gene_sorts[[gene]]
    tp <- 1 / (1 + log2(1 + tp))
    tp2 <- log2(1 + tp)
    tp_PPI_contribution[[gene]] <- tp2
    tp_PPI_contribution_norm[[gene]] <- tp2 / sum(tp2)
  }
  PPI_contribution[[PPI]] <- tp_PPI_contribution[results_of_PPIs[[PPI]]]
  PPI_contribution_norm[[PPI]] <- tp_PPI_contribution_norm[results_of_PPIs[[PPI]]]
}
rm(PPI)
all_considered_genes <- sort(unique(all_considered_genes))
gene_sort <- list()
for (PPI in PPIs)
{
  tsort <- (1:length(results_of_PPIs[[PPI]]))
  tsort <- (tsort / length(results_of_PPIs[[PPI]])) * log2(1 + c(1:length(results_of_PPIs[[PPI]])))
  dataFrame <- data.frame(gene = results_of_PPIs[[PPI]], sort = tsort)
  rownames(dataFrame) <- dataFrame$gene
  for (gene in results_of_PPIs[[PPI]])
  {
    tp <- dataFrame[gene, 2]
    names(tp) <- PPI
    gene_sort[[gene]] <- c(gene_sort[[gene]], tp)
  }
}
rm(PPI, dataFrame, gene, results_of_PPIs)

final_sort <- c()
for (gene in all_considered_genes)
{
  score <- mean(gene_sort[[gene]]) / length(gene_sort[[gene]])
  final_sort <- c(final_sort, score)
}
names(final_sort) <- all_considered_genes
final_sort <- sort(final_sort)
final_symbol <- names(final_sort)

contribute_ <- list()
contribute_norm <- list()
for (gene in all_considered_genes)
{
  tp <- gene_sort[[gene]]
  tp <- 1 / (1 + log2(1 + tp))
  tp2 <- log2(1 + tp)
  contribute_[[gene]] <- tp2
  contribute_norm[[gene]] <- tp2 / sum(tp2)
}
contribute_ <- contribute_[final_symbol]
contribute_norm <- contribute_norm[final_symbol]
saveRDS(contribute_norm, paste0("Results/Data/", cancer, "/Explaination_PPI_to_gene.rds"))
