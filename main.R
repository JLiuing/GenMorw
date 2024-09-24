## Prediction of driver gene prioritization for cancer patients and cohorts
# This script is responsible for identifying patient and cancer cohort driver Genes based on a Multi-omics random 
# walk with restart (GenMorw) model. In this model, the standard input is the GDC TCGA data provided by XenaBrowser 
# (https://xenabrowser.net/datapages/). Single Nucleotide Variations (SNVs) are required inputs, while others—including 
# Copy Number Variations (CNVs), DNA methylation data, gene expression data, and miRNA expression data—are optional. 
# Additionally, GenMorw comes with 12 PPI networks stored in the "Data/network_rds" directory, all of which are used 
# by default. Users can manually remove PPI files (.rds) from this directory or add their own networks for analysis.

# The following is the parameter statement:
# cores: Number of cores used for multi-threading computation.
# beta: Restart probability of the walking node in the RWR model.
# trans: Probability of the walking node transitioning between gene and patient layers in the heterogeneous network.
# guidance: guidance: Setting for the guidance pattern, where 0 indicates no guidance, 1 indicates standard guidance, and 2 indicates super guidance.
# GuidanceSet: Guidance set applicable when the guidance pattern is set to 1 or 2. Available options include: "cancermine",
#   "CancerSCEM", "CGC", "civic", "IntOGen", "MutPanning", "NCG", "OncoKB", "TSG" and "ALL".
# PPIs: The function AvailablePPI() automatically retrieves all PPIs from the "Data/network_rds" directory for constructing
#   the heterogeneous network. Users can freely combine PPIs, with options including "BioPlex", "CPDB", "HINT", "HumanNet", "InBioMap",
#   "IntAct", "IREF", "irefindex", "Mentha", "MULTINET", and "STRINGdb". For example, set "PPIs <- c("CPDB", "PCNET")". By default, all PPIs are used.
# retain_mmPPI: Logical value. True indicates that gene interaction information in the heterogeneous network will be retained, while False means 
#   it will not be retained. To construct the GenMorw-network, this must be set to True. Additionally, allocate 40-80 GB of disk space for each cancer.
# is_visualize: Logical value. True enables visualization of patient grouped plots for different data types, while False ignores visualization.

cores <- 4
beta <- 0.1
trans <- 0.01
guidance <- 1
GuidanceSet <- "NCG"
source("AvailablePPI.R")
PPIs <- AvailablePPI()
retain_mmPPI <- F
is_visualize <- T

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
