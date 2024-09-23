preprocess_fpkm <- function(path, cancer_name, is_scale_data = F, is_visualize = F) {
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  options(digits = 16)
  dir_temp <- paste0("Data/", cancer_name, "/fpkm")
  if (!dir.exists(dir_temp)) {
    dir.create(dir_temp)
  }

  fpkm_data <- fread(path, sep = "\t")
  tem <- fpkm_data[, 2]
  del_index <- tem %>% c()
  del_index <- del_index[[1]] %>% is.na()
  fpkm_data <- fpkm_data[!del_index, ] %>% na.omit()

  del_id <- which(duplicated(as.character(colnames(fpkm_data))))
  if (length(del_id) != 0) {
    fpkm_data <- fpkm_data[, -del_id]
  }

  thre <- (length(colnames(fpkm_data[, -1])) * 0.9) %>% ceiling()
  statis <- sapply(fpkm_data[, -1], function(x) {
    x == 0
  }) %>% rowSums()
  del_index <- which(statis >= thre)
  fpkm_data <- fpkm_data[-del_index, ]
  rm(del_index, thre, statis)

  ENSG_mapping <- read_tsv("Data/gencode.v22.annotation.gene.probeMap", show_col_types = F)
  ENSG_mapping <- ENSG_mapping[, c("id", "gene")]

  colnames(ENSG_mapping) <- c(colnames(fpkm_data)[1], colnames(ENSG_mapping)[2])
  fpkm_data <- inner_join(ENSG_mapping, fpkm_data, by = "Ensembl_ID")
  fpkm_data <- fpkm_data[, -1]

  fpkm_data <- fpkm_data %>%
    group_by(gene) %>%
    summarise(across(everything(), mean)) %>%
    as.data.frame()

  rownames(fpkm_data) <- fpkm_data$gene
  fpkm_data <- fpkm_data[, -1]
  rm(ENSG_mapping)

  cancer_sampleID <- colnames(fpkm_data)[!is.na(str_extract(colnames(fpkm_data), "^.+-0[1-9]A$"))]
  fpkm_data <- fpkm_data[, cancer_sampleID]
  nsample <- ncol(fpkm_data)
  rm(cancer_sampleID)

  suppressWarnings(fpkm_data_seu <- CreateSeuratObject(fpkm_data))
  var.num <- (length(rownames(fpkm_data)) * 0.1) %>% ceiling()
  var.num <- max(var.num, 2000)
  rm(fpkm_data)
  gc()

  fpkm_data_seu <- NormalizeData(fpkm_data_seu, normalization.method = "CLR", verbose = F)
  fpkm_data_seu <- FindVariableFeatures(fpkm_data_seu, nfeatures = var.num, verbose = F)
  var_genes <- VariableFeatures(fpkm_data_seu)
  write(var_genes, paste0(dir_temp, "/var_genes"))
  rm(var_genes)

  scale.genes <- rownames(fpkm_data_seu)
  fpkm_data_seu <- ScaleData(fpkm_data_seu, features = scale.genes, verbose = F)
  if (is_scale_data) {
    scale_data <- GetAssayData(fpkm_data_seu, slot = "scale.data", assay = "RNA") %>% as.data.frame()
    scale_data <- cbind(data.frame(gene = rownames(scale_data)), scale_data)
    fwrite(scale_data, paste0(dir_temp, "/scale_", cancer, ".csv.gz"))
    rm(scale_data)
  }

  if (nsample >= 60) {
    npcs <- 50 # pca default
  } else {
    npcs <- 30
  }

  fpkm_data_seu <- RunPCA(fpkm_data_seu, features = VariableFeatures(fpkm_data_seu), npcs = npcs, verbose = F)
  plot1 <- DimPlot(fpkm_data_seu, reduction = "pca", group.by = "orig.ident")
  plot2 <- ElbowPlot(fpkm_data_seu, ndims = (npcs - 10), reduction = "pca")
  plotc <- plot1 + plot2
  ggsave(paste0(dir_temp, "/pca.pdf"), plot = plotc, width = 8, height = 3.5)
  rm(plot1, plot2, plotc)

  pc.num <- 1:(npcs * 0.8)
  fpkm_data_seu <- FindNeighbors(fpkm_data_seu, dims = pc.num, prune.SNN = .01, k.param = (npcs * 0.8), n.trees = 100, verbose = F)
  fpkm_data_seu <- FindClusters(fpkm_data_seu, resolution = .95, algorithm = 2, verbose = F)
  metadata <- fpkm_data_seu@meta.data
  cell_cluster <- data.frame(patient_ID = rownames(metadata), cluster_ID = as.numeric(as.vector(metadata$seurat_clusters)))
  if (nsample >= 200) {
    sub_cell_cluster <- cell_cluster
    rownames(sub_cell_cluster) <- sub_cell_cluster$patient_ID
    sub_cell_cluster$cluster_ID <- 0
    groups <- unique(cell_cluster$cluster_ID) %>% sort()
    for (group in groups)
    {
      sub_fpkm_data_seu <- subset(fpkm_data_seu, cells = cell_cluster[cell_cluster$cluster_ID == group, 1])
      sub_fpkm_data_seu <- FindVariableFeatures(sub_fpkm_data_seu, nfeatures = var.num, verbose = F)
      scale.genes <- rownames(sub_fpkm_data_seu)
      sub_fpkm_data_seu <- ScaleData(sub_fpkm_data_seu, features = scale.genes, verbose = F)

      sub_npcs <- min(npcs, ncol(sub_fpkm_data_seu) - 1)
      sub_pc.num <- 1:floor(sub_npcs * 0.8)
      sub_fpkm_data_seu <- RunPCA(sub_fpkm_data_seu, features = VariableFeatures(sub_fpkm_data_seu), npcs = sub_npcs, verbose = F)
      sub_fpkm_data_seu <- FindNeighbors(sub_fpkm_data_seu, dims = sub_pc.num, prune.SNN = .01, k.param = floor(sub_npcs * 0.8), n.trees = 100, verbose = F)
      sub_fpkm_data_seu <- FindClusters(sub_fpkm_data_seu, resolution = .95, algorithm = 2, verbose = F)
      metadata <- sub_fpkm_data_seu@meta.data
      tcell_cluster <- data.frame(patient_ID = rownames(metadata), cluster_ID = as.numeric(as.vector(metadata$seurat_clusters)))
      rownames(tcell_cluster) <- tcell_cluster$patient_ID
      now_cluster_num <- sub_cell_cluster$cluster_ID %>%
        unique() %>%
        length()
      sub_cell_cluster[tcell_cluster$patient_ID, "cluster_ID"] <- tcell_cluster[tcell_cluster$patient_ID, "cluster_ID"] + now_cluster_num
    }
    rm(sub_fpkm_data_seu)
    sub_cell_cluster$cluster_ID <- sub_cell_cluster$cluster_ID - 1
    write.csv(sub_cell_cluster, paste0(dir_temp, "/patient_cluster.csv"), row.names = F)
  } else {
    write.csv(cell_cluster, paste0(dir_temp, "/patient_cluster.csv"), row.names = F)
  }
  SNN <- fpkm_data_seu@graphs$RNA_snn %>%
    as("matrix") %>%
    as.data.frame()
  fwrite(SNN, paste0(dir_temp, "/SNN.csv.gz"))
  rm(metadata, SNN)
  if (is_visualize) {
    try({
      fpkm_data_seu <- RunTSNE(fpkm_data_seu, dims = pc.num)
      plot1 <- DimPlot(fpkm_data_seu, reduction = "tsne")

      fpkm_data_seu <- RunUMAP(fpkm_data_seu, dims = pc.num, verbose = F)
      plot2 <- DimPlot(fpkm_data_seu, reduction = "umap")
      plotc <- plot1 + plot2
      ggsave(paste0(dir_temp, "/tsne_umap.pdf"), plot = plotc, width = 8, height = 3.5)
    })
  }
  rm(fpkm_data_seu)
}
