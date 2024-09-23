Sub_Model_Omics <- function(cancer, p,
                            beta = 0.1,
                            trans = 0.01,
                            scale_NN = 3,
                            base_NN = 5,
                            retain_mmPPI = F,
                            guidance = 1,
                            GuidanceSet = "NCG") {
  library(tidyverse)
  library(Matrix)
  options(digits = 16)
  gc()
  Sys.sleep(round(runif(1, 1, 30)))
  PPI_name <- p["PPI"]
  omic <- p["omic"]
  patient_clusters <- read.csv(paste0("Data/", cancer, "/", omic, "/patient_cluster.csv"))
  tpatients_correlation <- read_csv(paste0("Data/", cancer, "/", omic, "/SNN.csv.gz"), show_col_types = F) %>% as.data.frame()
  rownames(tpatients_correlation) <- colnames(tpatients_correlation)

  snv_cnv_path <- paste0("Data/", cancer, "/snv_combine_cnv.tsv.gz")
  if (file.exists(snv_cnv_path)) {
    patient_mut <- read_tsv(snv_cnv_path, show_col_types = F)
  } else {
    patient_mut <- read_tsv(paste0("Data/", cancer, "/mut_snv_filtered.tsv.gz"), show_col_types = F)
  }
  rm(snv_cnv_path)

  common_patients <- intersect(patient_clusters$patient_ID, colnames(patient_mut)[-1])
  clusters_ <- unique(patient_clusters$cluster_ID) %>% sort()
  patient_groups <- list()
  for (clust in clusters_)
  {
    nameID <- paste0("Group ", clust)
    input_patients <- intersect(patient_clusters[patient_clusters$cluster_ID == clust, 1], common_patients)
    if (length(input_patients) < 3) {
      next
    }
    patient_groups[[nameID]] <- input_patients
  }
  rm(clust, clusters_, patient_clusters, nameID, common_patients, input_patients)

  patient_mut_groups <- list()
  groupID <- names(patient_groups)
  for (group in groupID)
  {
    patientIDs <- patient_groups[[group]]
    mut_inf <- patient_mut[, c("gene", patientIDs)] %>% as.data.frame()
    rownames(mut_inf) <- mut_inf$gene
    mut_inf <- mut_inf[, -1]
    patient_mut_groups[[group]] <- as(as.matrix(mut_inf[as.logical(rowSums(mut_inf > 0)), ]), "sparseMatrix")
  }
  rm(group, patientIDs, patient_groups, mut_inf, patient_mut)

  PPI <- readRDS(paste0("Data/networks_rds/", PPI_name, ".rds"))
  col_names <- colnames(PPI)
  row_names <- rownames(PPI)
  PPI_genes <- intersect(col_names, row_names) %>% sort()
  PPI <- PPI[PPI_genes, PPI_genes]
  PPI <- PPI + t(PPI)
  PPI[PPI != 0] <- 1
  diag(PPI) <- 0
  rm(col_names, row_names, PPI_genes)
  PPI_groups <- list()
  All_PPI_genes <- rownames(PPI)
  for (group in groupID)
  {
    mut_inf <- patient_mut_groups[[group]]
    genes <- rownames(mut_inf)
    genes_ppi <- intersect(genes, All_PPI_genes)
    mut_inf <- mut_inf[genes_ppi, ]
    patient_mut_groups[[group]] <- mut_inf[, as.logical(colSums(mut_inf > 0))]
    PPI_groups[[group]] <- as(as.matrix(PPI[genes_ppi, genes_ppi]), "sparseMatrix")
  }
  rm(mut_inf, genes, genes_ppi, All_PPI_genes, PPI, group)

  patient_correlation_groups <- list()
  for (group in groupID)
  {
    patientIDs <- colnames(patient_mut_groups[[group]])
    patient_correlation_groups[[group]] <- tpatients_correlation[patientIDs, patientIDs] %>% as.matrix() - diag(length(patientIDs))
  }
  rm(group, tpatients_correlation, patientIDs)

  patient_correlation_groups_norm <- list()
  patient_mut_groups_norm <- list()
  t_patient_mut_groups_norm <- list()
  PPI_groups_norm <- list()

  for (group in groupID)
  {
    temp <- patient_correlation_groups[[group]] %>%
      apply(2, function(x) {
        if (sum(x)) {
          (x / sum(x)) * (1 - trans)
        } else {
          x
        }
      })

    patient_name_sort <- colnames(temp) %>% sort()
    patient_correlation_groups_norm[[group]] <- as(as.matrix(temp[patient_name_sort, patient_name_sort]), "sparseMatrix")

    temp <- patient_mut_groups[[group]] %>%
      apply(2, function(x) {
        if (sum(x)) {
          (x / sum(x)) * trans
        } else {
          x
        }
      })
    gene_name_sort <- rownames(temp) %>% sort()
    patient_mut_groups_norm[[group]] <- as(as.matrix(temp[gene_name_sort, patient_name_sort]), "sparseMatrix")

    temp <- t(patient_mut_groups[[group]]) %>%
      apply(2, function(x) {
        if (sum(x)) {
          (x / sum(x)) * trans
        } else {
          x
        }
      })
    t_patient_mut_groups_norm[[group]] <- as(as.matrix(temp[patient_name_sort, gene_name_sort]), "sparseMatrix")

    temp <- PPI_groups[[group]] %>%
      apply(2, function(x) {
        if (sum(x)) {
          (x / sum(x)) * (1 - trans)
        } else {
          x
        }
      })
    PPI_groups_norm[[group]] <- as(as.matrix(temp[gene_name_sort, gene_name_sort]), "sparseMatrix")
  }
  rm(
    group, PPI_groups, patient_correlation_groups,
    gene_name_sort, patient_name_sort, temp
  )

  prob_matrix <- list()
  for (group in groupID)
  {
    isolated_gene_data <- PPI_groups_norm[[group]] %>% colSums()
    isolated_genes <- isolated_gene_data[isolated_gene_data == 0] %>% names()

    temp1 <- cbind(PPI_groups_norm[[group]], patient_mut_groups_norm[[group]])
    temp2 <- cbind(t_patient_mut_groups_norm[[group]], patient_correlation_groups_norm[[group]])

    temp2[, isolated_genes] <- temp2[, isolated_genes] * (1 / trans)
    prob_matrix[[group]] <- as(as.matrix(rbind(temp1, temp2)), "sparseMatrix")
  }
  rm(
    group, temp1, temp2, isolated_gene_data, isolated_genes,
    patient_correlation_groups_norm, patient_mut_groups_norm, PPI_groups_norm, t_patient_mut_groups_norm
  )

  degs <- NULL
  if (file.exists(paste0("Data/DEGs_2022-04-27/", cancer, ".txt"))) {
    degs <- read_tsv(paste0("Data/DEGs_2022-04-27/", cancer, ".txt"), show_col_types = F)[, c("Gene Symbol", "Log2(Fold Change)")]
    degs$`Log2(Fold Change)` <- abs(degs$`Log2(Fold Change)`)
    degs <- degs %>%
      group_by(`Gene Symbol`) %>%
      summarise(`Log2(Fold Change)` = mean(`Log2(Fold Change)`)) %>%
      as.data.frame()
    rownames(degs) <- degs$`Gene Symbol`
  }

  influence_score <- list()
  sensibility_score <- list()
  source("Models/RWR.R")

  snv_cnv_path <- paste0("Data/", cancer, "/snv_combine_cnv.tsv.gz")
  if (file.exists(snv_cnv_path)) {
    patient_mut <- read_tsv(snv_cnv_path, show_col_types = F) %>% as.data.frame()
  } else {
    patient_mut <- read_tsv(paste0("Data/", cancer, "/mut_snv_filtered.tsv.gz"), show_col_types = F) %>% as.data.frame()
  }
  rownames(patient_mut) <- patient_mut$gene
  patient_mut$gene <- NULL
  patients_ <- colnames(patient_mut)
  rm(snv_cnv_path)
  gc()

  for (group in groupID)
  {
    ncol_ <- ncol(prob_matrix[[group]])
    data <- RWR(prob_matrix[[group]], beta, guidance = guidance, GuidanceSet = GuidanceSet)
    mut_data <- patient_mut_groups[[group]]
    patientIDs <- rownames(data)[!is.na(str_extract(rownames(data), "^TCGA.+$"))] %>% unique()
    gene_symbol <- setdiff(rownames(data), patientIDs)
    data_inf <- t(data[patientIDs, gene_symbol])
    data_sen <- data[gene_symbol, patientIDs]

    patients_sim <- data[patientIDs, patientIDs]
    if (retain_mmPPI) {
      path <- paste0("Models/Data/", cancer, "/mmPPI")
      if (!dir.exists(path)) {
        dir.create(path, recursive = T)
      }
      saveRDS(data[gene_symbol, gene_symbol], paste0(path, "/", PPI_name, "_", omic, "_", group, ".rds"))
    }
    rm(data)
    gc()
    diag(patients_sim) <- 0

    patient_inf <- list()
    patient_sen <- list()
    for (patient in patientIDs)
    {
      patient_score_inf <- list()
      patient_score_sen <- list()
      if (length(which(mut_data[, patient] != 0)) == 1) {
        patient_mut_genes <- which(mut_data[, patient] != 0) %>% names()
      } else {
        patient_mut_genes <- mut_data[which(mut_data[, patient] != 0), patient] %>%
          names() %>%
          sort()
      }

      tinf_score <- c()
      tsen_score <- c()

      tpatient_sim_inf <- patients_sim[patient, ] %>% sort(decreasing = T)
      tpatient_sim_sen <- patients_sim[, patient] %>% sort(decreasing = T)

      for (gene in patient_mut_genes)
      {
        if (length(degs) != 0) {
          if (gene %in% degs$`Gene Symbol`) {
            time_gene <- round(degs[gene, "Log2(Fold Change)"])

            ns1 <- base_NN + scale_NN * time_gene
            if (ns1 < length(tpatient_sim_inf)) {
              considered_patients_inf <- tpatient_sim_inf[1:ns1] %>% names()
              considered_patients_sen <- tpatient_sim_sen[1:ns1] %>% names()
            } else {
              considered_patients_inf <- tpatient_sim_inf %>% names()
              considered_patients_sen <- tpatient_sim_sen %>% names()
            }
          } else {
            if (base_NN < length(tpatient_sim_inf)) {
              considered_patients_inf <- tpatient_sim_inf[1:base_NN] %>% names()
              considered_patients_sen <- tpatient_sim_sen[1:base_NN] %>% names()
            } else {
              considered_patients_inf <- tpatient_sim_inf %>% names()
              considered_patients_sen <- tpatient_sim_sen %>% names()
            }
          }
        } else {
          if (base_NN < length(tpatient_sim_inf)) {
            considered_patients_inf <- tpatient_sim_inf[1:base_NN] %>% names()
            considered_patients_sen <- tpatient_sim_sen[1:base_NN] %>% names()
          } else {
            considered_patients_inf <- tpatient_sim_inf %>% names()
            considered_patients_sen <- tpatient_sim_sen %>% names()
          }
        }
        gene_mut <- patient_mut[gene, ] %>%
          t() %>%
          as.vector()
        names(gene_mut) <- patients_
        gene_mut <- gene_mut[gene_mut > 0]
        gene_patient <- names(gene_mut)
        considered_patients_inf <- intersect(considered_patients_inf, gene_patient)
        considered_patients_sen <- intersect(considered_patients_sen, gene_patient)

        if (length(considered_patients_inf) == 0 || length(considered_patients_sen) == 0) {
          gene_inf_score <- data_inf[gene, patient]
          gene_sen_score <- data_sen[gene, patient]
        } else {
          tt_sim_inf <- tpatient_sim_inf[considered_patients_inf]
          tt_sim_sen <- tpatient_sim_sen[considered_patients_sen]

          similar_weight_inf <- tt_sim_inf / sum(tt_sim_inf)
          similar_weight_sen <- tt_sim_sen / sum(tt_sim_sen)

          gene_inf_score <- data_inf[gene, patient] + sum(data_inf[gene, considered_patients_inf] * similar_weight_inf)
          gene_sen_score <- data_sen[gene, patient] + sum(data_sen[gene, considered_patients_sen] * similar_weight_sen)
        }
        tinf_score <- c(tinf_score, gene_inf_score)
        tsen_score <- c(tsen_score, gene_sen_score)
      }
      names(tinf_score) <- patient_mut_genes
      names(tsen_score) <- patient_mut_genes

      patient_score_inf[["mut_score"]] <- tinf_score
      patient_score_sen[["mut_score"]] <- tsen_score

      patient_inf[[patient]] <- patient_score_inf
      patient_sen[[patient]] <- patient_score_sen
    }
    influence_score[[group]] <- patient_inf
    sensibility_score[[group]] <- patient_sen
  }
  rm(
    group, data, mut_data, patientIDs, gene_symbol, data_inf, data_sen, patients_sim, patient_inf,
    patient_sen, patient, patient_score_inf, patient_score_sen, patient_mut_genes,
    patient_mut_groups, considered_patients_inf, considered_patients_sen,
    gene, gene_inf_score, gene_sen_score, similar_weight_inf, similar_weight_sen, time_gene,
    tinf_score, tpatient_sim_inf, tpatient_sim_sen, tsen_score, tt_sim_inf, tt_sim_sen
  )

  path <- paste0("Models/Data/", cancer, "/", PPI_name, "/", omic)
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  saveRDS(influence_score, paste0(path, "/influence_score.rds"))
  saveRDS(sensibility_score, paste0(path, "/sensibility_score.rds"))
}
