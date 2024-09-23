Cancer_name <- "GBM"

all_mut_genes <- scan(paste0("Data/", Cancer_name, "/all_mut_genes"), what = "c", quiet = T)
global_matrix <- matrix(rep(0, length(all_mut_genes) * length(all_mut_genes)), ncol = length(all_mut_genes))
colnames(global_matrix) <- all_mut_genes
rownames(global_matrix) <- all_mut_genes
mpath <- paste0("Models/Data/", Cancer_name, "/mmPPI")
all_group_files <- list.files(mpath)
i <- 0
for (file in all_group_files)
{
  i <- i + 1
  print(paste0(Cancer_name, " ", i, " / ", length(all_group_files)))
  demo <- readRDS(paste0(mpath, "/", file))
  diag(demo) <- 0
  demo <- apply(demo, 2, function(x) {
    thre <- quantile(x, 0.95)
    ifelse(x > thre, 1, 0)
  })

  global_matrix[rownames(demo), colnames(demo)] <- demo + global_matrix[rownames(demo), colnames(demo)]
}
global_matrix <- global_matrix[, colSums(global_matrix) != 0]
global_matrix <- global_matrix[rowSums(global_matrix) != 0, ]
saveRDS(global_matrix, paste0("Results/Data/", Cancer_name, "/", Cancer_name, "_GenMorw-network.rds"))
