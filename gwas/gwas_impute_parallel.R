#!/usr/bin/env Rscript
# gwas_impute_parallel.R — Parallel LFMM-based imputation using sNMF Q matrix
# For each SNP with missing genotypes, regresses observed values on Q (admixture
# proportions) and predicts missing values from the fitted model — equivalent to
# LEA's impute() but parallelised over SNP chunks.
#
# Usage: Rscript gwas_impute_parallel.R <lfmm_file> <run_dir> [n_cores] [chunk_size]
# Output: <lfmm_file without .lfmm>_imputed.lfmm  (and matching .012.pos symlink)

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

args       <- commandArgs(trailingOnly = TRUE)
f_lfmm     <- if (length(args) >= 1) args[1] else stop("Usage: gwas_impute_parallel.R <lfmm_file> <run_dir>")
run_dir    <- if (length(args) >= 2) args[2] else stop("Usage: gwas_impute_parallel.R <lfmm_file> <run_dir>")
n_cores    <- if (length(args) >= 3) as.integer(args[3]) else max(1L, detectCores() - 1L)
chunk_size <- if (length(args) >= 4) as.integer(args[4]) else 10000L

f_out    <- sub("\\.lfmm$", "_imputed.lfmm", f_lfmm)
f_pos_in <- sub("\\.lfmm$", ".012.pos", f_lfmm)
f_pos_out <- sub("\\.lfmm$", "_imputed.012.pos", f_lfmm)
q_rds    <- file.path(run_dir, "snmf_Q.rds")

cat("=== gwas_impute_parallel.R ===\n")
cat("Input   :", f_lfmm,  "\n")
cat("Q file  :", q_rds,   "\n")
cat("Output  :", f_out,   "\n")
cat("Cores   :", n_cores, "\n")
cat("Chunk   :", chunk_size, "SNPs\n\n")

# Load Q matrix (n_samples × K)
Q <- readRDS(q_rds)
cat("Q matrix:", nrow(Q), "x", ncol(Q), "\n")

# Read genotype matrix (samples × SNPs)
cat("Reading genotype matrix...\n")
Y <- as.matrix(fread(f_lfmm, header = FALSE))
cat("Matrix:", nrow(Y), "samples x", ncol(Y), "SNPs\n")
n_miss <- sum(Y == 9 | Y == -9)
cat("Missing genotypes:", n_miss,
    sprintf("(%.1f%% of all genotypes)\n\n", 100 * n_miss / length(Y)))

if (n_miss == 0L) {
  cat("No missing values — copying input to output.\n")
  file.copy(f_lfmm, f_out, overwrite = TRUE)
  file.copy(f_pos_in, f_pos_out, overwrite = TRUE)
  quit(save = "no")
}

# Recode missing as NA
Y_na <- Y
Y_na[Y == 9 | Y == -9] <- NA
rm(Y); gc()

# Parallelised imputation: one chunk of SNPs per worker
chunk_starts <- seq(1L, ncol(Y_na), by = chunk_size)
chunk_ends   <- pmin(chunk_starts + chunk_size - 1L, ncol(Y_na))
cat("Imputing", ncol(Y_na), "SNPs in", length(chunk_starts),
    "chunks on", n_cores, "cores...\n")

impute_chunk <- function(i) {
  idx <- chunk_starts[i]:chunk_ends[i]
  Yc  <- Y_na[, idx, drop = FALSE]

  for (j in seq_len(ncol(Yc))) {
    obs <- which(!is.na(Yc[, j]))
    mis <- which(is.na(Yc[, j]))
    if (length(mis) == 0L) next

    if (length(obs) < ncol(Q) + 1L) {
      # Too few observed for regression — fall back to observed mean
      Yc[mis, j] <- round(mean(Yc[obs, j], na.rm = TRUE))
      next
    }

    # OLS regression of observed genotypes on Q (with intercept)
    X_obs  <- cbind(1, Q[obs, , drop = FALSE])
    y_obs  <- Yc[obs, j]
    coef   <- tryCatch(
      solve(crossprod(X_obs), crossprod(X_obs, y_obs)),
      error = function(e) NULL
    )
    if (is.null(coef)) {
      Yc[mis, j] <- round(mean(y_obs))
      next
    }

    X_mis  <- cbind(1, Q[mis, , drop = FALSE])
    pred   <- as.numeric(X_mis %*% coef)
    Yc[mis, j] <- pmin(pmax(round(pred), 0L), 2L)
  }
  Yc
}

chunk_res <- mclapply(seq_along(chunk_starts), impute_chunk, mc.cores = n_cores)
Y_imp <- do.call(cbind, chunk_res)
rm(chunk_res); gc()

remaining <- sum(is.na(Y_imp))
cat("Remaining missing after imputation:", remaining, "\n")
if (remaining > 0L) {
  # Fallback: mean-impute any residual NAs
  col_means <- colMeans(Y_imp, na.rm = TRUE)
  for (j in which(colSums(is.na(Y_imp)) > 0L))
    Y_imp[is.na(Y_imp[, j]), j] <- round(col_means[j])
  cat("Residual NAs filled with column means.\n")
}

cat("Writing imputed matrix to:", f_out, "\n")
fwrite(as.data.table(Y_imp), f_out, sep = "\t", col.names = FALSE)

# Copy pos file so gwas_lfmm.R can find it
file.copy(f_pos_in, f_pos_out, overwrite = TRUE)
cat("Pos file copied to:", f_pos_out, "\n")
cat("Done.\n")
