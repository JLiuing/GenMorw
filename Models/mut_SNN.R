mut_SNN <- function(cancer, is_scale_data = F, is_visualize = F) {
  library(readr)
  library(tidyverse)
  library(Seurat)
  options(digits = 16)

  dir_temp <- paste0("Data/", cancer, "/mut")
  if (!dir.exists(dir_temp)) {
    dir.create(dir_temp)
  }

  snv_cnv_path <- paste0("Data/", cancer, "/snv_combine_cnv.tsv.gz")
  if (file.exists(snv_cnv_path)) {
    patient_mut <- read_tsv(snv_cnv_path, show_col_types = F) %>% as.data.frame()
  } else {
    patient_mut <- read_tsv(paste0("Data/", cancer, "/mut_snv_filtered.tsv.gz"), show_col_types = F) %>% as.data.frame()
  }

  if (length(rownames(patient_mut)) != length(rownames(na.omit(patient_mut)))) {
    patient_mut[is.na(patient_mut)] <- 0
  }

  rownames(patient_mut) <- patient_mut$gene
  patient_mut <- patient_mut[, -1]
  nsample <- ncol(patient_mut)

  suppressWarnings(mut_seu <- CreateSeuratObject(patient_mut))

  if (nsample >= 800) {
    var.num <- (length(rownames(patient_mut)) * 0.1) %>% ceiling()
  } else if (nsample < 800 && nsample >= 600) {
    var.num <- (length(rownames(patient_mut)) * 0.2) %>% ceiling()
  } else if (nsample < 600 && nsample >= 400) {
    var.num <- (length(rownames(patient_mut)) * 0.3) %>% ceiling()
  } else if (nsample < 400 && nsample >= 200) {
    var.num <- (length(rownames(patient_mut)) * 0.4) %>% ceiling()
  } else {
    var.num <- (length(rownames(patient_mut)) * 0.5) %>% ceiling()
  }

  rm(snv_cnv_path, patient_mut)
  gc()
  mut_seu <- NormalizeData(mut_seu, normalization.method = "CLR", verbose = F)
  mut_seu <- FindVariableFeatures(mut_seu, nfeatures = var.num, verbose = F)
  scale.genes <- rownames(mut_seu)
  mut_seu <- ScaleData(mut_seu, features = scale.genes, verbose = F)
  if (is_scale_data) {
    scale_data <- GetAssayData(mut_seu, slot = "scale.data", assay = "RNA") %>% as.data.frame()
    scale_data <- cbind(data.frame(gene = rownames(scale_data)), scale_data)
    write_csv(scale_data, paste0(dir_temp, "/scale_", cancer, ".csv.gz"))
    rm(scale_data)
  }
  rm(scale.genes)

  # PCA
  if (nsample >= 60) {
    npcs <- 50
  } else {
    npcs <- 30
  }
  mut_seu <- RunPCA(mut_seu, features = VariableFeatures(mut_seu), npcs = npcs, verbose = F)
  pc.num <- 1:(npcs * 0.6)
  mut_seu <- FindNeighbors(mut_seu, dims = pc.num, prune.SNN = .01, k.param = (npcs * 0.8), n.trees = 100, verbose = F)
  mut_seu <- FindClusters(mut_seu, resolution = .95, algorithm = 2, verbose = F)
  metadata <- mut_seu@meta.data
  cell_cluster <- data.frame(patient_ID = rownames(metadata), cluster_ID = as.numeric(as.vector(metadata$seurat_clusters)))
  if (nsample >= 200) {
    sub_cell_cluster <- cell_cluster
    rownames(sub_cell_cluster) <- sub_cell_cluster$patient_ID
    sub_cell_cluster$cluster_ID <- 0
    groups <- unique(cell_cluster$cluster_ID) %>% sort()
    for (group in groups)
    {
      sub_mut_seu <- subset(mut_seu, cells = cell_cluster[cell_cluster$cluster_ID == group, 1])
      sub_mut_seu <- FindVariableFeatures(sub_mut_seu, nfeatures = var.num, verbose = F)
      scale.genes <- rownames(sub_mut_seu)
      sub_mut_seu <- ScaleData(sub_mut_seu, features = scale.genes, verbose = F)

      sub_npcs <- min(npcs, ncol(sub_mut_seu) - 1)
      sub_pc.num <- 1:floor(sub_npcs * 0.8)
      sub_mut_seu <- RunPCA(sub_mut_seu, features = VariableFeatures(sub_mut_seu), npcs = npcs, verbose = F)
      sub_mut_seu <- FindNeighbors(sub_mut_seu, dims = sub_pc.num, prune.SNN = .01, k.param = floor(sub_npcs * 0.8), n.trees = 100, verbose = F)
      sub_mut_seu <- FindClusters(sub_mut_seu, resolution = .95, algorithm = 2, verbose = F)
      metadata <- sub_mut_seu@meta.data
      tcell_cluster <- data.frame(patient_ID = rownames(metadata), cluster_ID = as.numeric(as.vector(metadata$seurat_clusters)))
      rownames(tcell_cluster) <- tcell_cluster$patient_ID
      now_cluster_num <- sub_cell_cluster$cluster_ID %>%
        unique() %>%
        length()
      sub_cell_cluster[tcell_cluster$patient_ID, "cluster_ID"] <- tcell_cluster[tcell_cluster$patient_ID, "cluster_ID"] + now_cluster_num
    }
    rm(sub_mut_seu)
    sub_cell_cluster$cluster_ID <- sub_cell_cluster$cluster_ID - 1
    write.csv(sub_cell_cluster, paste0(dir_temp, "/patient_cluster.csv"), row.names = F)
  } else {
    write.csv(cell_cluster, paste0(dir_temp, "/patient_cluster.csv"), row.names = F)
  }
  SNN <- mut_seu@graphs$RNA_snn %>%
    as("matrix") %>%
    as.data.frame()
  write_csv(SNN, paste0(dir_temp, "/SNN.csv.gz"))
  rm(metadata, SNN)
  if (is_visualize) {
    try({
      mut_seu <- RunUMAP(mut_seu, dims = pc.num, verbose = F)
      plot <- DimPlot(mut_seu, reduction = "umap")
      ggsave(paste0(dir_temp, "/umap.pdf"), plot = plot, width = 6, height = 5.5)
    })
  }
  rm(mut_seu)
}
