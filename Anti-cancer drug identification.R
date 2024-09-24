# Estimation of potential anti-cancer drugs and database validation
#   Please execute this program after completing the prediction of driver genes for the cancer cohort. 
# The program is responsible for predicting the associations between the predicted driver genes and 
# drugs, and it validates these gene-drug pairs by searching across all databases on Entrez.
#   If you wish to validate the predicted genes for [cancer], ensure that the prioritization of 
# driver genes for the [cancer] cohort is completed, specifically confirming the existence of the file 
# "Results/Data/[cancer]/final_symbol_sort_patientFirst.txt". The associations and significance of 
# gene-drug pairs will be estimated and ultimately visualized at the end of the process.
# The following is the statement of the parameters:
#   Cancer_name: The TCGA abbreviation of the cancer to be analyzed
#   top_num: Number of top-ranked genes selected for gene-drug pair analysis
#   GDSC2: Drug sensitivity dataset (https://zenodo.org/records/5787145)
#   CCLE: Drug sensitivity dataset (https://zenodo.org/records/3905462)
#   cores: Number of cores used for multi-threading during searches on Entrez
#   retmax: Maximum number of search results returned from each database.

Cancer_name <- "LAML"
top_num <- 10
GDSC2 <- readRDS("path to GDSC2 dataset")
CCLE <- readRDS("path to CCLE dataset")
cores <- 4
retmax <- 20

library(rentrez)
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(PharmacoGx)
warnings(NULL)
predicted_genes <- read.table(paste0("Results/Data/", Cancer_name, "/final_symbol_sort_patientFirst.txt"))[1:top_num, 2]
GDSC2_drug_names <- drugNames(GDSC2)
CCLE_drug_names <- drugNames(CCLE)
GDSC2_drug_record <- data.frame()
CCLE_drug_record <- data.frame()
for (gene in predicted_genes)
{
  print(paste0(gene, " Analyzing."))
  feature_GDSC2 <- fNames(GDSC2, "rna")[which(featureInfo(GDSC2, "rna")$Symbol == gene)]
  if (length(feature_GDSC2) != 0) {
    sig.rna_GDSC2 <- drugSensitivitySig(
      object = GDSC2,
      mDataType = "rna",
      drugs = GDSC2_drug_names,
      features = feature_GDSC2,
      sensitivity.measure = "ic50_recomputed",
      molecular.summary.stat = "mean",
      sensitivity.summary.stat = "median",
      verbose = F
    )
    if (length(sig.rna_GDSC2) != 0) {
      call_data <- sig.rna_GDSC2@.Data[1, , ]
      new_rname <- paste0(gene, " + ", rownames(call_data))
      rownames(call_data) <- new_rname
      GDSC2_drug_record <- rbind(GDSC2_drug_record, call_data)
    }
  }
  feature_CCLE <- fNames(CCLE, "rna")[which(featureInfo(CCLE, "rna")$Symbol == gene)]
  if (length(feature_CCLE) != 0) {
    sig.rna_CCLE <- drugSensitivitySig(
      object = CCLE,
      mDataType = "rna",
      drugs = CCLE_drug_names,
      features = feature_CCLE,
      sensitivity.measure = "ic50_recomputed",
      molecular.summary.stat = "mean",
      sensitivity.summary.stat = "median",
      verbose = F
    )
    if (length(sig.rna_CCLE) != 0) {
      call_data <- sig.rna_CCLE@.Data[1, , ]
      new_rname <- paste0(gene, " + ", rownames(call_data))
      rownames(call_data) <- new_rname
      CCLE_drug_record <- rbind(CCLE_drug_record, call_data)
    }
  }
  print(paste0(gene, " Finished."))
}

GDSC2_drug_record <- GDSC2_drug_record[order(GDSC2_drug_record$fdr), ]
GDSC2_drug_record <- cbind(data.frame(Combination = rownames(GDSC2_drug_record)), GDSC2_drug_record)
CCLE_drug_record <- CCLE_drug_record[order(CCLE_drug_record$fdr), ]
CCLE_drug_record <- cbind(data.frame(Combination = rownames(CCLE_drug_record)), CCLE_drug_record)

all_dbs <- entrez_dbs()
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

myCluster <- makeCluster(cores)
registerDoParallel(myCluster)

dir_anti <- paste0("Results/Data/", Cancer_name, "/Anti-cancer drug records")
if (!dir.exists(dir_anti)) {
  dir.create(dir_anti)
}

combination_GDSC2 <- GDSC2_drug_record$Combination
foreach(i = 1:length(combination_GDSC2)) %dopar% {
  library(rentrez)
  search_entrez <- function(gene, drug, all_dbs) {
    search_term <- paste(gene, "AND", drug)
    results <- list()

    for (db in all_dbs) {
      cat("Searching in:", db, "\n")
      result <- tryCatch(
        {
          entrez_search(db = db, term = search_term, retmax = retmax)
        },
        error = function(e) {
          cat("Error in database:", db, "\n")
          return(NULL)
        }
      )
      results[[db]] <- result$ids
    }
    return(results)
  }
  comb <- combination_GDSC2[i]
  tp <- strsplit(comb, split = " + ", fixed = T)
  gene <- tp[[1]][1]
  drug <- tp[[1]][2]
  search_result <- search_entrez(gene, drug, all_dbs)
  search_result <- Filter(function(x) length(x) > 0, search_result)

  formatted_result <- sapply(names(search_result), function(db) {
    ids <- paste(search_result[[db]], collapse = ", ")
    paste0(db, " (", ids, ")")
  }, simplify = FALSE)
  final_output <- paste(formatted_result, collapse = ", ")
  saveRDS(final_output, paste0(dir_anti, "/", comb, ".rds"))
}

combination_CCLE <- CCLE_drug_record$Combination
foreach(i = 1:length(combination_CCLE)) %dopar% {
  library(rentrez)
  search_entrez <- function(gene, drug, all_dbs) {
    search_term <- paste(gene, "AND", drug)
    results <- list()

    for (db in all_dbs) {
      cat("Searching in:", db, "\n")
      result <- tryCatch(
        {
          entrez_search(db = db, term = search_term, retmax = 20) 
        },
        error = function(e) {
          cat("Error in database:", db, "\n")
          return(NULL)
        }
      )
      results[[db]] <- result$ids
    }
    return(results)
  }
  comb <- combination_CCLE[i]
  tp <- strsplit(comb, split = " + ", fixed = T)
  gene <- tp[[1]][1]
  drug <- tp[[1]][2]
  search_result <- search_entrez(gene, drug, all_dbs)
  search_result <- Filter(function(x) length(x) > 0, search_result)

  formatted_result <- sapply(names(search_result), function(db) {
    ids <- paste(search_result[[db]], collapse = ", ")
    paste0(db, " (", ids, ")")
  }, simplify = FALSE)
  final_output <- paste(formatted_result, collapse = ", ")
  saveRDS(final_output, paste0(dir_anti, "/", comb, ".rds"))
}
stopCluster(myCluster)

combination_GDSC2 <- GDSC2_drug_record$Combination
GDSC2_drug_record$evidence <- ""
for (comp in combination_GDSC2)
{
  tp <- readRDS(paste0(dir_anti, "/", comp, ".rds"))
  GDSC2_drug_record[comp, "evidence"] <- tp
}

combination_CCLE <- CCLE_drug_record$Combination
CCLE_drug_record$evidence <- ""
for (comp in combination_CCLE)
{
  tp <- readRDS(paste0(dir_anti, "/", comp, ".rds"))
  CCLE_drug_record[comp, "evidence"] <- tp
}

write.xlsx(GDSC2_drug_record, paste0("Results/Data/", Cancer_name, "/GDSC2_drug_record_with_evidence.xlsx"))
write.xlsx(CCLE_drug_record, paste0("Results/Data/", Cancer_name, "/CCLE_drug_record_with_evidence.xlsx"))
file.remove(list.files(dir_anti))

GDSC2_plot_df <- GDSC2_drug_record[, c(1, 2, 9)]
colnames(GDSC2_plot_df) <- c("Name", "Estimate", "FDR")
GDSC2_plot_df$FDR <- -log10(GDSC2_plot_df$FDR)
GDSC2_plot_df$FDR[GDSC2_plot_df$FDR >= 30] <- 30
GDSC2_plot_df$anno <- ""
GDSC2_plot_df[1:10, "anno"] <- GDSC2_plot_df[1:10, "Name"]

CCLE_plot_df <- CCLE_drug_record[, c(1, 2, 9)]
colnames(CCLE_plot_df) <- c("Name", "Estimate", "FDR")
CCLE_plot_df$FDR <- -log10(CCLE_plot_df$FDR)
CCLE_plot_df$FDR[CCLE_plot_df$FDR >= 10] <- 10
CCLE_plot_df$anno <- ""
CCLE_plot_df[1:10, "anno"] <- CCLE_plot_df[1:10, "Name"]

p_GDSC2 <- ggplot(data = GDSC2_plot_df, aes(x = Estimate, y = FDR)) +
  geom_point(aes(x = Estimate, y = FDR, color = FDR, size = FDR), alpha = 0.8) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  geom_text_repel(aes(label = anno, color = FDR),
    size = 3,
    max.overlaps = 100,
    key_glyph = draw_key_point
  ) +
  theme_classic() +
  labs(
    x = "Estimate",
    y = "-Log10(FDR)"
  )
ggsave(paste0("Results/Data/", Cancer_name, "/Gene_drug_within_GDSC2.pdf"), p_GDSC2, width = 6, height = 6)

p_CCLE <- ggplot(data = CCLE_plot_df, aes(x = Estimate, y = FDR)) +
  geom_point(aes(x = Estimate, y = FDR, color = FDR, size = FDR), alpha = 0.8) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")) +
  geom_text_repel(aes(label = anno, color = FDR),
    size = 3,
    max.overlaps = 100,
    key_glyph = draw_key_point
  ) +
  theme_classic() +
  labs(
    x = "Estimate",
    y = "-Log10(FDR)"
  )
ggsave(paste0("Results/Data/", Cancer_name, "/Gene_drug_within_CCLE.pdf"), p_CCLE, width = 6, height = 6)
