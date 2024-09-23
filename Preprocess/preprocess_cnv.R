preprocess_cnv <- function(path, cancer_name) {
  library(readr)
  library(tidyverse)
  options(digits = 16)

  dir_temp <- paste0("Data/", cancer_name)
  if (!dir.exists(dir_temp)) {
    dir.create(dir_temp)
  }

  mut_data <- read_tsv(path, show_col_types = F,col_names = F) %>% na.omit()
  del_id <- which(duplicated(as.character(mut_data[1,])))
  if(length(del_id)!=0)
  {
    mut_data <- mut_data[,-del_id]
  }
  colnames(mut_data) <- mut_data[1,]
  mut_data <- mut_data[-1,]
  
  ENSG_mapping <- read_tsv("Data/gencode.v22.annotation.gene.probeMap", show_col_types = F)
  ENSG_mapping <- ENSG_mapping[, c("id", "gene")]
  colnames(ENSG_mapping) <- c("Gene Symbol", "gene")

  CNV_data <- inner_join(ENSG_mapping, mut_data, by = "Gene Symbol")
  rm(mut_data, ENSG_mapping)
  CNV_data <- CNV_data[, -1]

  CNV_data2 <- sapply(CNV_data[, -1], function(x) {
    abs(as.integer(x))
  }) %>% as.data.frame()
  CNV_data <- cbind(data.frame(gene = CNV_data$gene), CNV_data2)
  rm(CNV_data2)
  CNV_data <- CNV_data[-which(duplicated(CNV_data$gene)),] %>% as.data.frame()

  write_tsv(CNV_data, paste0(dir_temp, "/mut_cnv.tsv.gz"))
  rm(CNV_data)
}
