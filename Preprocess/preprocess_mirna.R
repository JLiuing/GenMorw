preprocess_mirna <- function(path, cancer_name, is_scale_data = F, is_visualize = F) {
  library(readr)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  options(digits = 16)
  dir_temp <- paste0("Data/", cancer_name, "/mirna")
  if (!dir.exists(dir_temp)) {
    dir.create(dir_temp)
  }
  mirna_data <- read_tsv(path, show_col_types = F, col_names = F)

  del_id <- which(duplicated(as.character(mirna_data[1, ])))
  if (length(del_id) != 0) {
    mirna_data <- mirna_data[, -del_id]
  }
  colnames(mirna_data) <- mirna_data[1, ]
  mirna_data <- mirna_data[-1, ]
  MI <- mirna_data$miRNA_ID
  mirna_data <- sapply(mirna_data[, -1], function(x) {
    as.numeric(x)
  })
  mirna_data <- cbind(data.frame(`miRNA_ID` = MI), mirna_data)

  thre <- (length(colnames(mirna_data[, -1])) * 0.9) %>% ceiling()
  statis <- sapply(mirna_data[, -1], function(x) {
    x == 0
  }) %>% rowSums()
  del_index <- which(statis >= thre)
  mirna_data <- mirna_data[-del_index, ] %>% as.data.frame()
  rm(del_index, thre, statis)
  rownames(mirna_data) <- mirna_data[, 1]
  mirna_data <- mirna_data[, -1]

  cancer_sampleID <- colnames(mirna_data)[!is.na(str_extract(colnames(mirna_data), "^.+-0[1-9]A$"))]
  mirna_data <- mirna_data[, cancer_sampleID]
  nsample <- ncol(mirna_data)
  rm(cancer_sampleID, MI)

  mirna_data_seu <- CreateSeuratObject(mirna_data)
  var.num <- (length(rownames(mirna_data)) * 0.5) %>% ceiling()
  rm(mirna_data)
  gc()

  mirna_data_seu <- NormalizeData(mirna_data_seu, normalization.method = "CLR", verbose = F)
  mirna_data_seu <- FindVariableFeatures(mirna_data_seu, nfeatures = var.num, verbose = F)
  var_rnas <- VariableFeatures(mirna_data_seu)
  write(var_rnas, paste0(dir_temp, "/var_rnas"))
  rm(var_rnas)

  scale.rnas <- rownames(mirna_data_seu)
  mirna_data_seu <- ScaleData(mirna_data_seu, features = scale.rnas, verbose = F)

  if (is_scale_data) {
    scale_data <- GetAssayData(mirna_data_seu, slot = "scale.data", assay = "RNA") %>% as.data.frame()
    scale_data <- cbind(data.frame(mirna_ID = rownames(scale_data)), scale_data)
    write_csv(scale_data, paste0(dir_temp, "/scale_BRCA.csv.gz"))
    rm(scale_data)
  }

  if (nsample >= 60) {
    npcs <- 50
  } else {
    npcs <- 30
  }

  mirna_data_seu <- RunPCA(mirna_data_seu, features = VariableFeatures(mirna_data_seu), npcs = npcs, verbose = F)
  plot1 <- DimPlot(mirna_data_seu, reduction = "pca", group.by = "orig.ident")
  plot2 <- ElbowPlot(mirna_data_seu, ndims = (npcs - 10), reduction = "pca")
  plotc <- plot1 + plot2
  ggsave(paste0(dir_temp, "/pca.pdf"), plot = plotc, width = 8, height = 3.5)
  rm(plotc, plot1, plot2)

  pc.num <- 1:(npcs * 0.8)
  mirna_data_seu <- FindNeighbors(mirna_data_seu, dims = pc.num, prune.SNN = .01, k.param = (npcs * 0.8), n.trees = 100, verbose = F)
  mirna_data_seu <- FindClusters(mirna_data_seu, resolution = .95, algorithm = 2, verbose = F)
  metadata <- mirna_data_seu@meta.data
  cell_cluster <- data.frame(patient_ID = rownames(metadata), cluster_ID = as.numeric(as.vector(metadata$seurat_clusters)))
  if (nsample >= 200) {
    sub_cell_cluster <- cell_cluster
    rownames(sub_cell_cluster) <- sub_cell_cluster$patient_ID
    sub_cell_cluster$cluster_ID <- 0
    groups <- unique(cell_cluster$cluster_ID) %>% sort()
    for (group in groups)
    {
      sub_mirna_data_seu <- subset(mirna_data_seu, cells = cell_cluster[cell_cluster$cluster_ID == group, 1])
      sub_mirna_data_seu <- FindVariableFeatures(sub_mirna_data_seu, nfeatures = var.num, verbose = F)
      scale.genes <- rownames(sub_mirna_data_seu)
      sub_mirna_data_seu <- ScaleData(sub_mirna_data_seu, features = scale.genes, verbose = F)

      sub_npcs <- min(npcs, ncol(sub_mirna_data_seu) - 1)
      sub_pc.num <- 1:floor(sub_npcs * 0.8)
      sub_mirna_data_seu <- RunPCA(sub_mirna_data_seu, features = VariableFeatures(sub_mirna_data_seu), npcs = sub_npcs, verbose = F)
      sub_mirna_data_seu <- FindNeighbors(sub_mirna_data_seu, dims = sub_pc.num, prune.SNN = .01, k.param = floor(sub_npcs * 0.8), n.trees = 100, verbose = F)
      sub_mirna_data_seu <- FindClusters(sub_mirna_data_seu, resolution = .95, algorithm = 2, verbose = F)
      metadata <- sub_mirna_data_seu@meta.data
      tcell_cluster <- data.frame(patient_ID = rownames(metadata), cluster_ID = as.numeric(as.vector(metadata$seurat_clusters)))
      rownames(tcell_cluster) <- tcell_cluster$patient_ID
      now_cluster_num <- sub_cell_cluster$cluster_ID %>%
        unique() %>%
        length()
      sub_cell_cluster[tcell_cluster$patient_ID, "cluster_ID"] <- tcell_cluster[tcell_cluster$patient_ID, "cluster_ID"] + now_cluster_num
    }
    rm(sub_mirna_data_seu)
    sub_cell_cluster$cluster_ID <- sub_cell_cluster$cluster_ID - 1
    write.csv(sub_cell_cluster, paste0(dir_temp, "/patient_cluster.csv"), row.names = F)
  } else {
    write.csv(cell_cluster, paste0(dir_temp, "/patient_cluster.csv"), row.names = F)
  }
  SNN <- mirna_data_seu@graphs$RNA_snn %>%
    as("matrix") %>%
    as.data.frame()
  write_csv(SNN, paste0(dir_temp, "/SNN.csv.gz"))
  rm(metadata, SNN)

  if (is_visualize) {
    try({
      mirna_data_seu <- RunTSNE(mirna_data_seu, dims = pc.num)
      plot1 <- DimPlot(mirna_data_seu, reduction = "tsne")

      mirna_data_seu <- RunUMAP(mirna_data_seu, dims = pc.num, verbose = F)
      plot2 <- DimPlot(mirna_data_seu, reduction = "umap")
      plotc <- plot1 + plot2
      ggsave(paste0(dir_temp, "/tsne_umap.pdf"), plot = plotc, width = 8, height = 3.5)
    })
  }
  rm(mirna_data_seu)
}
