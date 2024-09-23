Predictions <- function(cancer, PPI) {
  library(tidyverse)
  options(digits = 16)

  mpath <- paste0("Models/Data/", cancer, "/", PPI)
  savepath <- paste0("Results/Data/", cancer, "/", PPI)
  if (!dir.exists(savepath)) {
    dir.create(savepath, recursive = T)
  }

  Harmonic_Mean <- function(data) {
    n_omic <- length(data)
    if (n_omic != 1) {
      gene_symbol <- names(data[[1]])
      gene_num <- length(gene_symbol)
      sum_ <- rep(0, gene_num)
      for (t in 1:n_omic)
      {
        sum_ <- sum_ + (data[[t]][gene_symbol])^-1
      }
      return(n_omic / sum_)
    } else {
      return(data[[1]])
    }
  }

  omics <- list.files(mpath)
  pred_inf <- list()
  pred_sen <- list()
  for (omic in omics)
  {
    path_omic_inf <- paste0(mpath, "/", omic, "/influence_score.rds")
    path_omic_sen <- paste0(mpath, "/", omic, "/sensibility_score.rds")

    # 取消分组
    temp_inf <- readRDS(path_omic_inf)
    tt_inf <- list()
    for (group in temp_inf)
    {
      tt_inf <- c(tt_inf, group)
    }
    pred_inf[[omic]] <- tt_inf

    temp_sen <- readRDS(path_omic_sen)
    tt_sen <- list()
    for (group in temp_sen)
    {
      tt_sen <- c(tt_sen, group)
    }
    pred_sen[[omic]] <- tt_sen
  }
  rm(omic, path_omic_inf, path_omic_sen, temp_inf, temp_sen, tt_inf, tt_sen, group)

  # 处理排序信息
  tpred_inf <- list()
  tpred_sen <- list()
  for (omic in omics)
  {
    # inf
    tomic_inf <- pred_inf[[omic]]
    store_inf <- list()
    patientIDs <- names(tomic_inf)
    for (patient in patientIDs)
    {
      temp <- (tomic_inf[[patient]]$mut_score)^-1
      store_inf[[patient]] <- (temp / sum(temp)) %>% sort(decreasing = T)
    }
    tpred_inf[[omic]] <- store_inf

    # sen
    tomic_sen <- pred_sen[[omic]]
    store_sen <- list()
    patientIDs <- names(tomic_sen)
    for (patient in patientIDs)
    {
      temp <- (tomic_sen[[patient]]$mut_score)
      store_sen[[patient]] <- (temp / sum(temp)) %>% sort(decreasing = T)
    }
    tpred_sen[[omic]] <- store_sen
  }
  pred_inf <- tpred_inf
  pred_sen <- tpred_sen
  rm(store_inf, store_sen, tomic_inf, tomic_sen, tpred_inf, tpred_sen, omic, patient, patientIDs, temp)
  saveRDS(pred_inf, paste0(savepath, "/pred_inf.rds"))
  saveRDS(pred_sen, paste0(savepath, "/pred_sen.rds"))

  inf_sen_integrated <- list()
  for (omic in omics)
  {
    omic_list <- list()
    omic_inf <- pred_inf[[omic]]
    omic_sen <- pred_sen[[omic]]
    # ID一致，取inf的ID
    patientIDs <- names(omic_inf)
    for (patient in patientIDs)
    {
      genes_patient <- names(omic_inf[[patient]])
      a <- omic_inf[[patient]][genes_patient]
      b <- omic_sen[[patient]][genes_patient]
      omic_list[[patient]] <- ((2 * a * b) / (a + b)) %>% sort(decreasing = T)
    }
    inf_sen_integrated[[omic]] <- omic_list
  }
  rm(omic_inf, omic_list, omic_sen, a, b, genes_patient, omic, patient, patientIDs)
  saveRDS(inf_sen_integrated, paste0(savepath, "/inf_sen_integrated.rds"))

  # 将信息汇总到病人上
  patients_omics <- list()
  patientIDs <- names(inf_sen_integrated$mut)
  for (patient in patientIDs)
  {
    patient_omic <- list()
    for (omic in omics)
    {
      if (patient %in% names(inf_sen_integrated[[omic]])) {
        temp <- inf_sen_integrated[[omic]][[patient]]
        patient_omic[[omic]] <- temp / sum(temp)
      }
    }
    patients_omics[[patient]] <- patient_omic
  }
  rm(patient_omic, patient, omic)

  all_considered_genes <- c()
  patient_results <- list()
  patient_results_symbol <- list()
  for (patient in patientIDs)
  {
    patient_results[[patient]] <- patients_omics[[patient]] %>%
      Harmonic_Mean() %>%
      sort(decreasing = T)
    patient_results_symbol[[patient]] <- names(patient_results[[patient]])
    all_considered_genes <- c(all_considered_genes, names(patient_results[[patient]]))
  }
  all_considered_genes <- all_considered_genes %>%
    unique() %>%
    sort()
  rm(patient, Harmonic_Mean)
  saveRDS(patient_results, paste0(savepath, "/patient_results.rds"))
  saveRDS(patient_results_symbol, paste0(savepath, "/patient_results_symbol.rds"))

  gene_sorts <- list()
  for (gene in all_considered_genes)
  {
    gene_sorts[[gene]] <- c()
  }
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
  }
  rm(gene, patient, dataFrame, mutgene)

  all_considered_genes <- all_considered_genes %>% sort()
  final_sort <- c()
  for (gene in all_considered_genes)
  {
    score <- mean(gene_sorts[[gene]]) / length(gene_sorts[[gene]])
    final_sort <- c(final_sort, score)
  }
  rm(gene, score)
  names(final_sort) <- all_considered_genes
  final_sort <- final_sort %>% sort()

  saveRDS(final_sort, paste0(savepath, "/cohort_sort.rds"))
  final_symbol <- final_sort %>% names()

  write.table(final_symbol, paste0(savepath, "/cohort_symbol.txt"), quote = F, col.names = F)
}
