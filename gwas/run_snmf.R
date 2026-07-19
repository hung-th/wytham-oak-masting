#!/usr/bin/env Rscript
# run_snmf.R — Run sNMF on LD-pruned genotype matrix and save Q matrix
# Usage: Rscript run_snmf.R <run_dir> <lfmm_file> [n_cores] [k_max]
# Example: Rscript run_snmf.R 1-vcf/original 1-vcf/original/WWO_37WGS.lfmm 100 1

suppressPackageStartupMessages({ library(LEA); library(data.table) })

args      <- commandArgs(trailingOnly = TRUE)
run_dir   <- if (length(args) >= 1) args[1] else stop("Usage: run_snmf.R <run_dir> <lfmm_file>")
lfmm_in   <- if (length(args) >= 2) args[2] else stop("Usage: run_snmf.R <run_dir> <lfmm_file>")
n_cores   <- if (length(args) >= 3) as.integer(args[3]) else max(1L, parallel::detectCores() - 1L)
k_max     <- if (length(args) >= 4) as.integer(args[4]) else 1L

prefix    <- gsub("\\.lfmm$", "", basename(lfmm_in))
snmf_lfmm <- file.path(run_dir, paste0(prefix, "_snmf_sub.lfmm"))
out_rds   <- file.path(run_dir, "snmf_Q.rds")

cat("=== run_snmf.R ===\n")
cat("Run dir  :", run_dir,   "\n")
cat("Input    :", lfmm_in,   "\n")
cat("Prefix   :", prefix,    "\n")
cat("Cores    :", n_cores,   "\n")
cat("K max    :", k_max,     "\n\n")

snmf_proj_path <- file.path(run_dir, paste0(prefix, "_snmf_sub.snmfProject"))

if (file.exists(snmf_proj_path)) {
  cat("Reloading existing sNMF project...\n")
  proj <- load.snmfProject(snmf_proj_path)
} else {
  cat("Reading lfmm for sNMF subsample...\n")
  mat <- as.matrix(fread(lfmm_in, header = FALSE))
  cat("Matrix:", nrow(mat), "samples x", ncol(mat), "SNPs\n")
  set.seed(225)
  sub_idx <- sample(ncol(mat), min(100000L, ncol(mat)))
  fwrite(as.data.table(mat[, sub_idx]), snmf_lfmm, sep = "\t", col.names = FALSE)
  rm(mat); gc()

  cat("Running sNMF K=1..", k_max, "(3 reps)...\n")
  proj <- snmf(snmf_lfmm, K = seq_len(k_max), entropy = TRUE,
               CPU = n_cores, repetitions = 3, project = "new")
}

mean_ce <- sapply(seq_len(k_max), function(k)
  mean(sapply(1:3, function(r)
    tryCatch(cross.entropy(proj, K = k, run = r), error = function(e) NA_real_))))
cat("Cross-entropy by K:", paste(round(mean_ce, 4), collapse = "  "), "\n")

# Save Q for K=1 best run (used by gwas_lfmm.R)
ce_k1    <- sapply(1:3, function(r)
  tryCatch(cross.entropy(proj, K = 1, run = r), error = function(e) NA_real_))
best_run <- which.min(ce_k1)
q_f <- file.path(run_dir, paste0(prefix, "_snmf_sub.snmf"),
                 "K1", paste0("run", best_run),
                 paste0(prefix, "_snmf_sub_r", best_run, ".1.Q"))
Q_mat <- as.matrix(read.table(q_f, header = FALSE))
saveRDS(Q_mat, out_rds)
cat("Q matrix (K=1, run", best_run, "):", nrow(Q_mat), "x", ncol(Q_mat), "\n")
cat("Saved:", out_rds, "\n")
