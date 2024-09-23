# main program
GenMorw <- function(cancer, omics_list, PPIs,
                    cores = 6,
                    beta = 0.1,
                    trans = 0.01,
                    scale_NN = 3,
                    base_NN = 5,
                    retain_mmPPI = F,
                    guidance = 1,
                    is_visualize = F,
                    GuidanceSet = "NCG") {
  options(digits = 16)
  Sys.setenv("VROOM_CONNECTION_SIZE" = 99999999)
  library(tidyverse)
  if (file.exists(omics_list$mut_snv)) {
    source("Preprocess/preprocess_snv.R")
    print("Preprocessing snv data.")
    preprocess_snv(path = omics_list$mut_snv, cancer_name = cancer) %>%
      suppressMessages() %>%
      suppressPackageStartupMessages()
    gc()

    if (file.exists(omics_list$mut_cnv)) {
      source("Preprocess/preprocess_cnv.R")
      print("Preprocessing cnv data.")
      preprocess_cnv(path = omics_list$mut_cnv, cancer_name = cancer) %>%
        suppressMessages() %>%
        suppressPackageStartupMessages()
      gc()

      source("Preprocess/snv_combine_cnv.R")
      print("Combining cnv to snv.")
      snv_combine_cnv(cancer_name = cancer) %>%
        suppressMessages() %>%
        suppressPackageStartupMessages()
      gc()
    } else {
      print("There is no cnv data!")
    }
  } else {
    print("Please input SNV mutation data!")
    return()
  }

  if (file.exists(omics_list$methy)) {
    print("Preprocessing methy data and grouping patients.")
    source("Preprocess/preprocess_methy.R")
    preprocess_methy(path = omics_list$methy, cancer_name = cancer, is_visualize = is_visualize) %>%
      suppressMessages() %>%
      suppressPackageStartupMessages()
    gc()
  } else {
    print("There is no methy data!")
  }

  if (file.exists(omics_list$fpkm)) {
    source("Preprocess/preprocess_fpkm.R")
    print("Preprocessing fpkm data and grouping patients.")
    preprocess_fpkm(path = omics_list$fpkm, cancer_name = cancer, is_visualize = is_visualize) %>%
      suppressMessages() %>%
      suppressPackageStartupMessages()
    gc()
  } else {
    print("There is no fpkm data!")
  }

  if (file.exists(omics_list$mirna)) {
    source("Preprocess/preprocess_mirna.R")
    print("Preprocessing mirna data and grouping patients.")
    preprocess_mirna(path = omics_list$mirna, cancer_name = cancer, is_visualize = is_visualize) %>%
      suppressMessages() %>%
      suppressPackageStartupMessages()
    gc()
  } else {
    print("There is no mirna data!")
  }

  dir_temp <- paste0("Models/Data")
  if (!dir.exists(dir_temp)) {
    dir.create(dir_temp)
  }
  dir_temp2 <- paste0("Results/Data")
  if (!dir.exists(dir_temp)) {
    dir.create(dir_temp2)
  }
  rm(dir_temp, dir_temp2)

  source("Models/mut_SNN.R")
  print("Grouping patients based on mutation data.")
  mut_SNN(cancer = cancer, is_visualize = is_visualize)
  gc()

  omics <- list.dirs(paste0("Data/", cancer), full.names = F, recursive = F)

  suppressPackageStartupMessages(library(foreach))
  suppressPackageStartupMessages(library(doParallel))

  myCluster <- makeCluster(cores)
  registerDoParallel(myCluster)
  nPPI <- length(PPIs)
  nomic <- length(omics)
  para <- list()
  tindex <- 0
  for (i in 1:nPPI)
  {
    for (j in 1:nomic)
    {
      tindex <- tindex + 1
      temp <- c(PPIs[i], omics[j])
      names(temp) <- c("PPI", "omic")
      para[[tindex]] <- temp
    }
  }
  rm(i, j, tindex, temp)

  print("Performing RWR on heterogeneous networks.")
  nn <- nPPI * nomic
  MES <- foreach(p = 1:nn) %dopar% {
    source("Models/Sub_Model_Omics.R")
    Sub_Model_Omics(
      cancer = cancer,
      p = para[[p]],
      beta = beta,
      trans = trans,
      scale_NN = scale_NN,
      base_NN = base_NN,
      retain_mmPPI = retain_mmPPI,
      guidance = guidance,
      GuidanceSet = GuidanceSet
    )
    mes <- "OK!"
  }

  print("Integrating results.")
  MES_pred <- foreach(pi = 1:nPPI) %dopar% {
    source("Results/Predictions.R")
    Predictions(cancer = cancer, PPI = PPIs[pi])
    mes <- "OK!"
  }
  stopCluster(myCluster)

  results_of_PPIs <- list()
  all_considered_genes <- c()
  Results_mpath <- paste0("Results/Data")
  for (PPI in PPIs)
  {
    rpath <- paste0(Results_mpath, "/", cancer, "/", PPI, "/cohort_symbol.txt")
    results_of_PPIs[[PPI]] <- read.table(rpath)[, 2]
    all_considered_genes <- c(all_considered_genes, results_of_PPIs[[PPI]])
  }
  rm(PPI, rpath)
  all_considered_genes <- sort(unique(all_considered_genes))
  gene_sort <- list()
  for (gene in all_considered_genes)
  {
    gene_sort[[gene]] <- c()
  }

  for (PPI in results_of_PPIs)
  {
    tsort <- (1:length(PPI))
    tsort <- (tsort / length(PPI)) * log2(1 + c(1:length(PPI)))
    dataFrame <- data.frame(gene = PPI, sort = tsort)
    rownames(dataFrame) <- dataFrame$gene
    for (gene in PPI)
    {
      gene_sort[[gene]] <- c(gene_sort[[gene]], dataFrame[gene, 2])
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

  write.table(final_symbol, paste0(Results_mpath, "/", cancer, "/final_symbol_sort_patientFirst.txt"), quote = F, col.names = F)
  rm(gene, score, gene_sort, final_sort, final_symbol)

  patient_sort_symbol <- list()
  for (PPI in PPIs)
  {
    tpatients <- readRDS(paste0(Results_mpath, "/", cancer, "/", PPI, "/patient_results_symbol.rds"))
    patientsID <- names(tpatients)
    for (patient in patientsID)
    {
      pdata <- tpatients[[patient]]
      if (length(pdata) == 0) {
        next
      }
      tsort <- 1:length(pdata)
      tsort <- (tsort / length(pdata)) * log2(1 + c(1:length(pdata)))
      sort_frame <- data.frame(gene = pdata, sort = tsort)
      rownames(sort_frame) <- sort_frame$gene
      if (patient %in% names(patient_sort_symbol)) {
        for (gene in pdata)
        {
          if (gene %in% names(patient_sort_symbol[[patient]])) {
            patient_sort_symbol[[patient]][[gene]] <- c(patient_sort_symbol[[patient]][[gene]], sort_frame[gene, 2])
          } else {
            patient_sort_symbol[[patient]][[gene]] <- c(sort_frame[gene, 2])
          }
        }
      } else {
        patient_sort_symbol[[patient]] <- list()
        for (gene in pdata)
        {
          patient_sort_symbol[[patient]][[gene]] <- c(sort_frame[gene, 2])
        }
      }
    }
  }
  rm(PPI, tpatients, patientsID, patient, pdata, sort_frame, gene)

  patients_symbol <- list()
  patientID <- names(patient_sort_symbol)
  all_genes <- c()
  for (patient in patientID)
  {
    tdata <- patient_sort_symbol[[patient]]
    pgenes <- names(tdata)
    all_genes <- c(all_genes, pgenes)
    pgenes_score <- c()
    for (gene in pgenes)
    {
      score <- mean(tdata[[gene]]) / (length(tdata[[gene]]))
      pgenes_score <- c(pgenes_score, score)
    }
    names(pgenes_score) <- pgenes
    pgenes_score <- pgenes_score %>% sort()
    patients_symbol[[patient]] <- pgenes_score %>% names()
  }
  rm(patient, tdata, pgenes, pgenes_score, gene, score, patient_sort_symbol)
  saveRDS(patients_symbol, paste0(Results_mpath, "/", cancer, "/single_patients_pred.rds"))
}
