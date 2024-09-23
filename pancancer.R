cancers <- c(
  "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH",
  "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG",
  "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
)

gene_sort <- list()
for(cancer in cancers)
{
  path <- paste0("Results/Data/",cancer,"/single_patients_pred.rds")
  cohort <- readRDS(path)
  for(patient in cohort)
  {
    tsort <- 1:length(patient)
    tsort <- (tsort / length(patient)) * log2(1 + c(1:length(patient)))
    dataFrame <- data.frame(gene=patient,sort=tsort)
    rownames(dataFrame) <- dataFrame$gene
    for(gene in patient)
    {
      gene_sort[[gene]] <- c(gene_sort[[gene]],dataFrame[gene,"sort"])
    }
  }
}

all_considered_genes <- sort(names(gene_sort))
final_sort <- c()
for(gene in all_considered_genes)
{
  score <- mean(gene_sort[[gene]])/length(gene_sort[[gene]])
  final_sort <- c(final_sort,score)
}
names(final_sort) <- all_considered_genes
final_sort <- sort(final_sort)
final_symbol <- names(final_sort)
write(final_symbol,"Results/Data/pancancer_prediction.txt")
