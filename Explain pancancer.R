## Retrieval of contribution percentages of various cancer cohorts to pancancer prioritized genes
# After completing the identification of driver genes across multiple cancer cohorts, GenMorw can 
# integrate these cohorts into a pan-cancer prioritization. The program can analyze the contribution
# composition of each gene in the prioritization, with the contributions derived collectively from
# multiple cancers. By default, the program uses all cancer data available under the "Results/Data" 
# directory for analysis. Users can also modify the "cancers" parameter to customize settings.

cancers <- list.files("Results/Data")

gene_sort <- list()
for (cancer in cancers)
{
  path <- paste0("Results/Data/", cancer, "/single_patients_pred.rds")
  cohort <- readRDS(path)
  for (patient in cohort)
  {
    tsort <- 1:length(patient)
    tsort <- (tsort / length(patient)) * log2(1 + c(1:length(patient)))
    dataFrame <- data.frame(gene = patient, sort = tsort)
    rownames(dataFrame) <- dataFrame$gene
    for (gene in patient)
    {
      tp <- c(dataFrame[gene, "sort"])
      names(tp) <- cancer
      gene_sort[[gene]] <- c(gene_sort[[gene]], tp)
    }
  }
}

all_considered_genes <- sort(names(gene_sort))
final_sort <- c()
for (gene in all_considered_genes)
{
  score <- mean(gene_sort[[gene]]) / length(gene_sort[[gene]])
  final_sort <- c(final_sort, score)
}
names(final_sort) <- all_considered_genes
final_sort <- sort(final_sort)
final_symbol <- names(final_sort)

cancer_contribution <- list()
for (gene in all_considered_genes)
{
  tp <- gene_sort[[gene]]
  tp <- 1 / (1 + log2(1 + tp))
  tp2 <- aggregate(tp, by = list(cancer = names(tp)), FUN = sum)
  tp3 <- log2(1 + tp2$x)
  names(tp3) <- tp2$cancer
  tp3 <- tp3 / sum(tp3)
  cancer_contribution[[gene]] <- tp3
}
cancer_contribution <- cancer_contribution[final_symbol]
saveRDS(cancer_contribution, "Results/Data/Cancer_contribution.rds")
