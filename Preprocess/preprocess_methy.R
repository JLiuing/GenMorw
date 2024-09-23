preprocess_methy <- function(path, cancer_name, is_scale_data = F, is_visualize = F) {
  library(data.table)
  library(tidyverse)
  library(Seurat)
  library(Matrix)
  dir_temp <- paste0("Data/", cancer_name, "/methy")
  if (!dir.exists(dir_temp)) {
    dir.create(dir_temp)
  }

  methy_data <- fread(path, sep = "\t")
  tem <- methy_data[, 2]
  del_index <- tem %>% c()
  del_index <- del_index[[1]] %>% is.na()
  methy_data <- methy_data[!del_index, ] %>% na.omit()
  rm(tem, del_index)

  del_id <- which(duplicated(as.character(colnames(methy_data))))
  if (length(del_id) != 0) {
    methy_data <- methy_data[, -del_id]
  }

  mapping <- read_table("Data/illuminaMethyl450_hg38_GDC", show_col_types = F)
  mapping <- mapping[mapping$gene != ".", ]
  mapping <- mapping[!grepl(",", mapping$gene), ] %>% as.data.frame()
  rownames(mapping) <- mapping$`#id`
  mapping <- mapping[, c(1, 2)]

  colnames(mapping) <- c(colnames(methy_data)[1], colnames(mapping)[2])
  methy_data <- inner_join(mapping, methy_data, by = "Composite Element REF")
  methy_data <- methy_data[, -1]

  methy_data <- methy_data %>%
    group_by(gene) %>%
    summarise(across(everything(), mean)) %>%
    as.data.frame()

  rownames(methy_data) <- methy_data$gene
  methy_data <- methy_data[, -1]

  cancer_sampleID <- colnames(methy_data)[!is.na(str_extract(colnames(methy_data), "^.+-0[1-9]A$"))]
  methy_data <- methy_data[, cancer_sampleID]
  nsample <- ncol(methy_data)
  rm(cancer_sampleID)

  suppressWarnings(methy_data_seu <- CreateSeuratObject(methy_data))
  var.num <- (length(rownames(methy_data)) * 0.1) %>% ceiling()
  var.num <- max(var.num, 2000)
  rm(methy_data, mapping)
  gc()

  methy_data_seu <- NormalizeData(methy_data_seu, normalization.method = "CLR", verbose = F)
  methy_data_seu <- FindVariableFeatures(methy_data_seu, nfeatures = var.num, verbose = F)
  var_genes <- VariableFeatures(methy_data_seu)
  write(var_genes, paste0(dir_temp, "/var_genes"))
  rm(var_genes)

  scale.genes <- rownames(methy_data_seu)
  methy_data_seu <- ScaleData(methy_data_seu, features = scale.genes, verbose = F)
  if (is_scale_data) {
    scale_data <- GetAssayData(methy_data_seu, slot = "scale.data", assay = "RNA") %>% as.data.frame()
    scale_data <- cbind(data.frame(gene = rownames(scale_data)), scale_data)
    fwrite(scale_data, paste0(dir_temp, "/scale_", cancer, ".csv.gz"))
    rm(scale_data)
  }

  if (nsample >= 60) {
    npcs <- 50
  } else {
    npcs <- 30
  }

  methy_data_seu <- RunPCA(methy_data_seu, features = VariableFeatures(methy_data_seu), npcs = npcs, verbose = F)
  plot1 <- DimPlot(methy_data_seu, reduction = "pca", group.by = "orig.ident")
  plot2 <- ElbowPlot(methy_data_seu, ndims = (npcs - 10), reduction = "pca")
  plotc <- plot1 + plot2
  ggsave(paste0(dir_temp, "/pca.pdf"), plot = plotc, width = 8, height = 3.5)

  pc.num <- 1:(npcs * 0.8)
  methy_data_seu <- FindNeighbors(methy_data_seu, dims = pc.num, prune.SNN = .01, k.param = (npcs * 0.8), n.trees = 100, verbose = F)
  methy_data_seu <- FindClusters(methy_data_seu, resolution = .95, algorithm = 2, verbose = F)
  metadata <- methy_data_seu@meta.data
  cell_cluster <- data.frame(patient_ID = rownames(metadata), cluster_ID = as.numeric(as.vector(metadata$seurat_clusters)))

  if (nsample >= 200) {
    sub_cell_cluster <- cell_cluster
    rownames(sub_cell_cluster) <- sub_cell_cluster$patient_ID
    sub_cell_cluster$cluster_ID <- 0
    groups <- unique(cell_cluster$cluster_ID) %>% sort()
    for (group in groups)
    {
      sub_methy_data_seu <- subset(methy_data_seu, cells = cell_cluster[cell_cluster$cluster_ID == group, 1])
      sub_methy_data_seu <- FindVariableFeatures(sub_methy_data_seu, nfeatures = var.num, verbose = F)
      scale.genes <- rownames(sub_methy_data_seu)
      sub_methy_data_seu <- ScaleData(sub_methy_data_seu, features = scale.genes, verbose = F)

      sub_npcs <- min(npcs, ncol(sub_methy_data_seu) - 1)
      sub_pc.num <- 1:floor(sub_npcs * 0.8)
      sub_methy_data_seu <- RunPCA(sub_methy_data_seu, features = VariableFeatures(sub_methy_data_seu), npcs = sub_npcs, verbose = F)
      sub_methy_data_seu <- FindNeighbors(sub_methy_data_seu, dims = sub_pc.num, prune.SNN = .01, k.param = floor(sub_npcs * 0.8), n.trees = 100, verbose = F)
      sub_methy_data_seu <- FindClusters(sub_methy_data_seu, resolution = .95, algorithm = 2, verbose = F)
      metadata <- sub_methy_data_seu@meta.data
      tcell_cluster <- data.frame(patient_ID = rownames(metadata), cluster_ID = as.numeric(as.vector(metadata$seurat_clusters)))
      rownames(tcell_cluster) <- tcell_cluster$patient_ID
      now_cluster_num <- sub_cell_cluster$cluster_ID %>%
        unique() %>%
        length()
      sub_cell_cluster[tcell_cluster$patient_ID, "cluster_ID"] <- tcell_cluster[tcell_cluster$patient_ID, "cluster_ID"] + now_cluster_num
    }
    rm(sub_methy_data_seu)
    sub_cell_cluster$cluster_ID <- sub_cell_cluster$cluster_ID - 1
    write.csv(sub_cell_cluster, paste0(dir_temp, "/patient_cluster.csv"), row.names = F)
  } else {
    write.csv(cell_cluster, paste0(dir_temp, "/patient_cluster.csv"), row.names = F)
  }

  SNN <- methy_data_seu@graphs$RNA_snn %>%
    as("matrix") %>%
    as.data.frame()
  fwrite(SNN, paste0(dir_temp, "/SNN.csv.gz"))
  rm(metadata, SNN)

  if (is_visualize) {
    try({
      methy_data_seu <- RunTSNE(methy_data_seu, dims = pc.num)
      plot1 <- DimPlot(methy_data_seu, reduction = "tsne")

      methy_data_seu <- RunUMAP(methy_data_seu, dims = pc.num, verbose = F)
      plot2 <- DimPlot(methy_data_seu, reduction = "umap")
      plotc <- plot1 + plot2
      ggsave(paste0(dir_temp, "/tsne_umap.pdf"), plot = plotc, width = 8, height = 3.5)
    })
  }
  rm(methy_data_seu)
}
