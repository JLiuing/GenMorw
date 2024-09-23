snv_combine_cnv <- function(cancer_name) {
  library(tidyverse)
  library(readr)
  options(digits = 16)

  path_snv <- paste0("Data/", cancer_name, "/mut_snv_filtered.tsv.gz")
  path_cnv <- paste0("Data/", cancer_name, "/mut_cnv.tsv.gz")
  SNV_data <- read_tsv(path_snv,show_col_types = F) %>% as.data.frame()
  CNV_data <- read_tsv(path_cnv,show_col_types = F) %>% as.data.frame()
  rownames(SNV_data) <- SNV_data$gene
  rownames(CNV_data) <- CNV_data$gene

  aids <- colnames(CNV_data)
  ids <- aids[!is.na(str_extract(aids, "^TCGA.+\\..+$"))]
  if (length(ids) != 0) {
    ids_split_after <- str_sub(ids, 1, 16)
    t <- 0
    for (i in ids)
    {
      t <- t + 1
      aids[which(aids == i)] <- ids_split_after[t]
    }
    colnames(CNV_data) <- aids
    rm(i, t, ids, ids_split_after, aids)
  }

  patients <- colnames(SNV_data)[-1]
  patients_cnv <- colnames(CNV_data)[-1]

  for (p in patients)
  {
    snv_p <- SNV_data[, c("gene", p)]
    snv_p <- filter(snv_p, snv_p[, 2] != 0)
    rownames(snv_p) <- snv_p[, 1]

    if (p %in% patients_cnv) {
      cnv_p <- CNV_data[, c("gene", p)]
      cnv_p <- filter(cnv_p, cnv_p[, 2] != 0)
      cnvs <- sum(cnv_p[, p])

      if (cnvs != 0) {
        snv_cnv_p <- intersect(snv_p[, 1], cnv_p[, 1])
        if (length(snv_cnv_p) != 0) {
          ratio <- 2 / (1 + exp(-log10(1 + cnvs)))

          replace_data <- snv_p[snv_cnv_p, 2] * ratio
          replace_data <- lapply(replace_data, function(x) {
            min(1, x)
          }) %>% unlist()
          names(replace_data) <- snv_cnv_p

          SNV_data[snv_cnv_p, p] <- replace_data[snv_cnv_p]
        }
      }
    }
  }
  rm(cnv_p, snv_p, snv_cnv_p, p, ratio, replace_data, cnvs)
  write_tsv(SNV_data, paste0("Data/", cancer_name, "/snv_combine_cnv.tsv.gz"))
  rm(SNV_data)
}
