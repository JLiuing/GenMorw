cores <- 4 # Set the cores
beta <- 0.1
trans <- 0.01
guidance <- 1 # 0: 无引导   1:引导    2:超级引导
GuidanceSet <- "NCG"
source("AvailablePPI.R")
PPIs <- AvailablePPI()
retain_mmPPI <- F
is_visualize <- T

# Considered cancers
cancers <- list.files("Cancer Data")
cancers <- unlist(lapply(strsplit(cancers, split = "GDC_"), function(x) {
  x[[2]]
}))
source("GenMorw.R")
for (cancer in cancers)
{
  omics_list <- list(
    mut_snv = paste0("Cancer Data/GDC_", cancer, "/TCGA-", cancer, ".mutect2_snv.tsv.gz"),
    mut_cnv = paste0("Cancer Data/GDC_", cancer, "/TCGA-", cancer, ".gistic.tsv.gz"),
    methy = paste0("Cancer Data/GDC_", cancer, "/TCGA-", cancer, ".methylation450.tsv.gz"),
    fpkm = paste0("Cancer Data/GDC_", cancer, "/TCGA-", cancer, ".htseq_fpkm.tsv.gz"),
    mirna = paste0("Cancer Data/GDC_", cancer, "/TCGA-", cancer, ".mirna.tsv.gz")
  )

  GenMorw(
    cancer = cancer,
    omics_list = omics_list,
    PPIs = PPIs,
    cores = cores,
    beta = beta,
    trans = trans,
    guidance = guidance,
    GuidanceSet = GuidanceSet,
    retain_mmPPI = retain_mmPPI,
    is_visualize=is_visualize
  )
}
