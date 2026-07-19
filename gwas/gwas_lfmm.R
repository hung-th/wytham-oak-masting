#!/usr/bin/env Rscript
# gwas_lfmm.R  —  parallel LFMM2 GWAS
# Usage: Rscript gwas_lfmm.R <run_dir> <data_csv> [n_cores] [chunk_size] [k_gwas] [excl_sample] [input_lfmm]
# Example (n37 full):   Rscript gwas_lfmm.R 1-vcf/n37 phenotypes/WWO_37WGS_data.csv 32 500000 1
# Example (n36):        Rscript gwas_lfmm.R 1-vcf/n36 phenotypes/WWO_37WGS_data.csv 32 500000 1 Ox002104
# Example (n37 pruned): Rscript gwas_lfmm.R 1-vcf/n37 phenotypes/WWO_37WGS_data.csv 32 500000 1 "" 1-vcf/n37/WWO_37WGS_pruned.lfmm
#
# input_lfmm (arg 7): optional path to override the auto-detected imputed .lfmm file.
#   Use this to run on the LD-pruned SNP set instead of the full imputed set.
#   The matching .012.pos file must exist at the same path with .lfmm → .012.pos substitution.
#   Outputs are written with suffix "_pruned" to avoid overwriting full-set results:
#   lfmm2_fit_pruned.rds, lfmm2_res_pruned.rds, lfmm.pruned.{p,q,z}.tsv
#
# Outputs (written to <run_dir>/):
#   lfmm2_fit[_pruned].rds  — fitted LFMM2 model (latent factors)
#   lfmm2_res[_pruned].rds  — combined p, z, gif, q per trait (list)
#   lfmm[.pruned].{p,q,z}.tsv

suppressPackageStartupMessages({
  library(LEA)
  library(data.table)
  library(parallel)
  library(qvalue)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args       <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript gwas_lfmm.R <run_dir> <data_csv> [n_cores] [chunk_size] [k_gwas] [excl_sample] [input_lfmm]")

run_dir      <- args[1]
data_csv     <- args[2]
n_cores      <- if (length(args) >= 3) as.integer(args[3]) else max(1L, detectCores() - 1L)
chunk_size   <- if (length(args) >= 4) as.integer(args[4]) else 500000L
k_gwas       <- if (length(args) >= 5) as.integer(args[5]) else 1L
excl_sample  <- if (length(args) >= 6) args[6] else ""
input_lfmm   <- if (length(args) >= 7 && nzchar(args[7])) args[7] else ""

# ── File detection ─────────────────────────────────────────────────────────────
lfmm_files <- list.files(run_dir, pattern = "\\.lfmm$")
lfmm_base  <- lfmm_files[!grepl("pruned|imputed|snmf", lfmm_files)]
PREFIX     <- gsub("\\.lfmm$", "", lfmm_base[1])

pruned_mode <- nzchar(input_lfmm)

if (pruned_mode) {
  f_imputed <- input_lfmm
  f_pos     <- sub("\\.lfmm$", ".012.pos", input_lfmm)
  # Derive suffix from input filename so different pruned sets get separate checkpoints
  sfx_base  <- gsub("^.*?([^/\\\\]+)\\.lfmm$", "\\1", input_lfmm)
  sfx_base  <- sub(paste0("^", PREFIX, "_?"), "", sfx_base)
  if (!nzchar(sfx_base)) sfx_base <- "pruned"
  ck_fit    <- file.path(run_dir, paste0("lfmm2_fit_", sfx_base, ".rds"))
  ck_res    <- file.path(run_dir, paste0("lfmm2_res_", sfx_base, ".rds"))
  tsv_sfx   <- sfx_base
} else {
  f_imputed <- file.path(run_dir, paste0(PREFIX, "_imputed.lfmm"))
  f_pos     <- file.path(run_dir, paste0(PREFIX, ".012.pos"))
  ck_fit    <- file.path(run_dir, "lfmm2_fit.rds")
  ck_res    <- file.path(run_dir, "lfmm2_res.rds")
  tsv_sfx   <- ""
}

snmf_dir   <- file.path(run_dir, paste0(PREFIX, "_snmf_sub.snmf"))
snmf_pfx   <- paste0(PREFIX, "_snmf_sub")

cat("=== gwas_lfmm.R ===\n")
cat("Run dir    :", run_dir,    "\n")
cat("Prefix     :", PREFIX,     "\n")
cat("Mode       :", if (pruned_mode) "LD-pruned" else "full imputed", "\n")
cat("Input      :", f_imputed,  "\n")
cat("Pos file   :", f_pos,      "\n")
cat("Data CSV   :", data_csv,   "\n")
cat("Cores      :", n_cores,    "\n")
cat("Chunk size :", chunk_size, "\n")
cat("K (LFMM2) :", k_gwas,     "\n")
cat("Excl sample:", if (nzchar(excl_sample)) excl_sample else "(none)", "\n\n")

stopifnot(
  "input lfmm not found" = file.exists(f_imputed),
  "pos file not found"   = file.exists(f_pos),
  "data csv not found"   = file.exists(data_csv)
)

# ── Phenotype matrix X ─────────────────────────────────────────────────────────
cat("[1/4] Preparing phenotype matrix...\n")
raw <- read.csv(data_csv, header = TRUE)

if (nzchar(excl_sample)) {
  raw <- raw[raw$GENOME != excl_sample, ]
  cat("Excluded sample:", excl_sample, "| Remaining:", nrow(raw), "rows\n")
}

# Soil texture PCA (sand/silt/clay are compositional)
soil_pca     <- prcomp(raw[, c("SAND", "SILT", "CLAY")], center = TRUE, scale. = TRUE)
raw$SOIL_PC1 <- soil_pca$x[, 1]

GWAS_VARS <- c("VIZ_COUNT", "MATURE_ACORNS", "IMMAT_ACORNS", "ENLARGED_CUPS",
               "FLOWERS", "CANOPY_CLOSURE", "SPRING_PHENO", "SOIL_PC1", "MIDNOV_LAI")

raw[c("MATURE_ACORNS", "IMMAT_ACORNS")] <- log1p(raw[c("MATURE_ACORNS", "IMMAT_ACORNS")])

X <- as.matrix(raw[, GWAS_VARS])
cat("X matrix:", nrow(X), "samples x", ncol(X), "traits\n\n")

# ── Q matrix from sNMF binary files ───────────────────────────────────────────
cat("[2/4] Loading Q matrix from sNMF...\n")
ck_Q <- file.path(run_dir, "snmf_Q.rds")

if (file.exists(ck_Q)) {
  Q_mat <- readRDS(ck_Q)
  cat("Loaded Q from checkpoint:", ck_Q, "\n")
} else if (dir.exists(snmf_dir)) {
  # LEA writes .Q files as space-separated text (not binary) — use read.table
  read_Q_txt <- function(path)
    as.matrix(read.table(path, header = FALSE))
  read_ce_bin <- function(path)
    tryCatch(slot(dget(path), "crossEntropy"), error = function(e) NA_real_)

  ce_rows <- list()
  Q_runs  <- list()
  for (r in 1:3) {
    q_f  <- file.path(snmf_dir, paste0("K", k_gwas), paste0("run", r),
                      paste0(snmf_pfx, "_r", r, ".", k_gwas, ".Q"))
    sc_f <- file.path(snmf_dir, paste0("K", k_gwas), paste0("run", r),
                      paste0(snmf_pfx, "_r", r, ".", k_gwas, ".snmfClass"))
    if (!file.exists(q_f)) next
    ce_rows[[r]] <- data.frame(rep = r, ce = read_ce_bin(sc_f))
    Q_runs[[r]]  <- read_Q_txt(q_f)
  }
  ce_tbl   <- do.call(rbind, ce_rows)
  best_run <- ce_tbl$rep[which.min(ce_tbl$ce)]
  Q_mat    <- Q_runs[[best_run]]
  cat("Best run at K =", k_gwas, ": run", best_run,
      "(CE =", round(min(ce_tbl$ce, na.rm=TRUE), 4), ")\n")
  saveRDS(Q_mat, ck_Q)
} else {
  stop("No sNMF project directory found at ", snmf_dir,
       "\nRun gwas_impute.R first, or provide ", ck_Q)
}
cat("Q matrix:", nrow(Q_mat), "x", ncol(Q_mat), "\n\n")

# ── Read imputed genotype matrix ───────────────────────────────────────────────
cat("[3/4] Reading imputed genotype matrix...\n")
Y    <- as.matrix(fread(f_imputed, header = FALSE))
loci <- fread(f_pos, header = FALSE, col.names = c("CHROM", "POS"))
n_samples <- nrow(Y)
n_snps    <- ncol(Y)
cat("Matrix:", n_samples, "samples x", n_snps, "SNPs\n\n")

stopifnot("Row mismatch: Y vs X" = nrow(Y) == nrow(X),
          "Row mismatch: Y vs Q" = nrow(Y) == nrow(Q_mat))

# ── Mean-impute any remaining missing values (coded 9) ────────────────────────
# The full imputed set has no 9s; the pruned set may have a small number.
# Replace 9s with column mean (2 × allele frequency) — standard for few missing.
missing_mask <- Y == 9 | Y == -9 | is.na(Y)
n_miss <- sum(missing_mask)
if (n_miss > 0L) {
  cat("Mean-imputing", n_miss, "missing genotypes in input matrix...\n")
  Y_na <- Y
  Y_na[missing_mask] <- NA
  col_means <- colMeans(Y_na, na.rm = TRUE)
  rm(Y_na)
  # Vectorised replacement — no per-column loop
  Y[missing_mask] <- round(col_means[col(Y)[missing_mask]])
  cat("Done.\n\n")
}

# ── LFMM2 fit ─────────────────────────────────────────────────────────────────
cat("[4/4] LFMM2...\n")

if (file.exists(ck_fit)) {
  cat("Loading fitted model from checkpoint:", ck_fit, "\n")
  lfmm2_fit <- readRDS(ck_fit)
} else {
  cat("Fitting LFMM2 (K =", k_gwas, ")...\n")
  lfmm2_fit <- lfmm2(input = Y, env = X, K = k_gwas)
  saveRDS(lfmm2_fit, ck_fit)
  cat("Fit saved to:", ck_fit, "\n")
}

# ── Parallel lfmm2.test by SNP chunk ─────────────────────────────────────────
if (file.exists(ck_res)) {
  cat("Loading results from checkpoint:", ck_res, "\n")
  res <- readRDS(ck_res)
} else {
  chunk_starts <- seq(1L, n_snps, by = chunk_size)
  chunk_ends   <- pmin(chunk_starts + chunk_size - 1L, n_snps)
  n_chunks     <- length(chunk_starts)
  cat("Testing", n_snps, "SNPs in", n_chunks, "chunks on", n_cores, "cores\n")
  cat("(genomic.control=FALSE per chunk; GIF applied globally after combining)\n\n")

  run_chunk <- function(i) {
    idx   <- chunk_starts[i]:chunk_ends[i]
    chunk <- Y[, idx, drop = FALSE]
    # genomic.control=FALSE — GIF must be computed globally
    lfmm2.test(lfmm2_fit, input = chunk, env = X, genomic.control = FALSE)
  }

  if (.Platform$OS.type == "unix") {
    chunk_res <- mclapply(seq_len(n_chunks), run_chunk, mc.cores = n_cores)
  } else {
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("lfmm2_fit", "X", "Y", "chunk_starts", "chunk_ends"),
                  envir = environment())
    clusterEvalQ(cl, library(LEA))
    chunk_res <- parLapply(cl, seq_len(n_chunks), run_chunk)
    stopCluster(cl)
  }

  # Combine chunks: zscores are traits × SNPs — cbind across chunks
  # Use Z-scores (not p-values) to recompute p-values globally, matching the
  # internal GIF computation in lfmm2.test(genomic.control=TRUE).
  z_raw <- do.call(cbind, lapply(chunk_res, `[[`, "zscores"))  # traits × SNPs
  rm(chunk_res, Y); gc()

  # Recompute p-values from Z-scores globally (two-sided chi-squared)
  chi2_raw <- z_raw^2  # traits × SNPs

  # ── Global genomic control per trait ─────────────────────────────────────
  cat("Applying global genomic control per trait...\n")
  chi2_null_median <- qchisq(0.5, df = 1, lower.tail = FALSE)
  n_traits <- nrow(z_raw)

  gif    <- numeric(n_traits)
  p_gc   <- chi2_raw
  q_mat  <- chi2_raw

  for (t in seq_len(n_traits)) {
    gif[t]      <- median(chi2_raw[t, ], na.rm = TRUE) / chi2_null_median
    chi2_corr   <- chi2_raw[t, ] / gif[t]
    p_gc[t, ]   <- pchisq(chi2_corr, df = 1, lower.tail = FALSE)
    q_mat[t, ]  <- qvalue(p_gc[t, ])$qvalues
  }

  gif_df <- data.frame(Trait = GWAS_VARS, GIF = round(gif, 4))
  cat("GIF per trait:\n"); print(gif_df); cat("\n")

  # Transpose to SNPs × traits for output
  p_out <- t(p_gc)
  q_out <- t(q_mat)
  z_out <- t(z_raw)
  colnames(p_out) <- colnames(q_out) <- colnames(z_out) <- GWAS_VARS

  res <- list(p = p_out, q = q_out, z = z_out, gif = gif_df, loci = loci)
  saveRDS(res, ck_res)
  cat("Results saved to:", ck_res, "\n")

  # Write flat TSV files (suffix distinguishes pruned from full-set results)
  sfx <- if (nzchar(tsv_sfx)) paste0(".", tsv_sfx) else ""
  write.table(cbind(loci, p_out), file.path(run_dir, paste0("lfmm", sfx, ".p.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(loci, q_out), file.path(run_dir, paste0("lfmm", sfx, ".q.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(loci, z_out), file.path(run_dir, paste0("lfmm", sfx, ".z.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(gif_df, file.path(run_dir, paste0("lfmm", sfx, "_gif.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

sig <- rowSums(res$q < 0.05) >= 1
cat("\n=== Summary ===\n")
cat("Total significant SNPs (Q<0.05 in >=1 trait):", sum(sig), "\n")
print(data.frame(Trait=GWAS_VARS, Sig_SNPs=colSums(res$q < 0.05)))
cat("\nDone. Download", ck_res, "to continue analysis locally.\n")
