RWR <- function(prob_matrix, beta, iter_time_max = 1000, gap = 1e-9, guidance = 1, GuidanceSet = "NCG") {
  library(Matrix)
  NGs <- scan("Data/NGs", what = "c", quiet = T)
  Benchmark <- scan(paste0("Benchmark Genes/", GuidanceSet), what = "c", quiet = T)

  ncol_ <- ncol(prob_matrix)
  iter0 <- diag(ncol_)
  colnames(iter0) <- colnames(prob_matrix)
  rownames(iter0) <- rownames(prob_matrix)

  if (guidance == 1) {
    NGs <- intersect(NGs, colnames(iter0))
    if (length(NGs) != 0) {
      diag(iter0[NGs, NGs]) <- round(2 / (1 + 1.01^(-length(NGs))), 6)
    }

    Benchmark <- intersect(Benchmark, colnames(iter0))
    if (length(Benchmark) != 0) {
      diag(iter0[Benchmark, Benchmark]) <- round((2 + log10(length(Benchmark))) / (2 * (1 + log10(length(Benchmark)))), 6)
    }
  } else if (guidance == 2) {
    NGs <- intersect(NGs, colnames(iter0))
    if (length(NGs) != 0) {
      diag(iter0[NGs, NGs]) <- 10
    }

    Benchmark <- intersect(Benchmark, colnames(iter0))
    if (length(Benchmark) != 0) {
      diag(iter0[Benchmark, Benchmark]) <- 0.1
    }
  }

  iter1 <- iter0
  iter0 <- as(iter0, "sparseMatrix")
  gc()
  for (i in 1:iter_time_max)
  {
    iter2 <- (1 - beta) * prob_matrix %*% iter1 + beta * iter0
    if (i %% 5 == 0) {
      mean_err <- sqrt(sum((iter2 - iter1)^2)) / (ncol_ * ncol_)
      if (mean_err <= gap) {
        return(iter2)
      }
    }
    iter1 <- iter2
  }
  return(iter2)
}
