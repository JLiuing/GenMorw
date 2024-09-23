preprocess_snv <- function(path, cancer_name, remain_synonymous_variant = T) {
  library(readr)
  library(tidyverse)
  library(Matrix)
  dir_temp <- paste0("Data/", cancer_name)
  if (!dir.exists(dir_temp)) {
    dir.create(dir_temp)
  }
  options(digits = 16)
  mut_data <- read_tsv(path, show_col_types = F)
  samples <- mut_data %>%
    group_by(Sample_ID) %>%
    summarise(n = n(), gene_num = length(gene))
  del_id <- samples$Sample_ID[which(samples$gene_num > 800)]
  mut_data_filtered <- filter(mut_data, !Sample_ID %in% del_id)
  rm(samples, del_id)
  if (!remain_synonymous_variant) {
    mut_data_filtered <- filter(mut_data_filtered, effect != "synonymous_variant")
    mut_data_filtered <- filter(mut_data_filtered, effect != "splice_region_variant;synonymous_variant")
    mut_data_filtered <- filter(mut_data_filtered, effect != "synonymous_variant;NMD_transcript_variant")
  }

  mut_matrix_filtered <- mut_data_filtered[, c("Sample_ID", "gene", "dna_vaf")]
  rm(mut_data_filtered)
  ngene_filtered <- mut_matrix_filtered$gene %>%
    unique() %>%
    length()
  nsample_filtered <- mut_matrix_filtered$Sample_ID %>%
    unique() %>%
    length()

  all_mut_genes <- mut_matrix_filtered$gene %>%
    unique()
  write(all_mut_genes, paste0(dir_temp, "/all_mut_genes"))

  mapping_samples_filtered <- data.frame(Index = 1:nsample_filtered, Sample_ID = sort(unique(mut_matrix_filtered$Sample_ID)))
  mapping_genes_filtered <- data.frame(Index = 1:ngene_filtered, Symbol = sort(unique(mut_matrix_filtered$gene)))
  rownames(mapping_genes_filtered) <- mapping_genes_filtered$Symbol
  rownames(mapping_samples_filtered) <- mapping_samples_filtered$Sample_ID

  mut_matrix_filtered <- mut_matrix_filtered %>%
    group_by(Sample_ID, gene) %>%
    summarise(mean_vaf = mean(dna_vaf)) %>%
    as.data.frame()
  mut_matrix_filtered <- cbind(mut_matrix_filtered,
    sample_index = mapping_samples_filtered[mut_matrix_filtered$Sample_ID, "Index"],
    gene_index = mapping_genes_filtered[mut_matrix_filtered$gene, "Index"]
  )
  mut_sparseMatirx_filtered <- sparseMatrix(mut_matrix_filtered$gene_index,
    mut_matrix_filtered$sample_index,
    x = mut_matrix_filtered$mean_vaf,
    dimnames = list(genes = mapping_genes_filtered$Symbol, samples = mapping_samples_filtered$Sample_ID)
  )
  rm(mut_matrix_filtered, ngene_filtered, nsample_filtered)
  dense_mut_filtered <- as(mut_sparseMatirx_filtered, "matrix") %>% as.data.frame()
  rm(mut_sparseMatirx_filtered)

  dense_mut_filtered <- cbind(data.frame(gene = rownames(dense_mut_filtered)), dense_mut_filtered)
  rm(mapping_genes_filtered, mapping_samples_filtered)
  write_tsv(dense_mut_filtered, paste0(dir_temp, "/mut_snv_filtered.tsv.gz"))
  rm(dense_mut_filtered)
}
