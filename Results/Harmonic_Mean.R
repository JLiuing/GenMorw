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
