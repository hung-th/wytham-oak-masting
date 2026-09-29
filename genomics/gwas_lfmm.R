#!/usr/bin/env Rscript
# gwas_lfmm.R — genome-wide association tests with LFMM2 (LEA)
#
# Usage:
#   Rscript genomics/gwas_lfmm.R <genotypes.lfmm> <positions.012.pos> <variables.csv> <out.rds> \
#                                [chunk_size] [K] [n_cores]
# As run for the study (from the project root):
#   Rscript genomics/gwas_lfmm.R \
#     data/raw/genomics/genotypes/WWO_37WGS_pruned_imputed.lfmm \
#     data/raw/genomics/genotypes/WWO_37WGS_pruned_imputed.012.pos \
#     data/derived/gwas_variables.csv cache/gwas/lfmm2_results.rds
#
# Inputs
#   genotypes.lfmm     imputed genotype matrix (samples x SNPs; 0/1/2), rows in the order of
#                      WWO_37WGS_pruned.012.indv (output of gwas_impute_parallel.R)
#   positions.012.pos  SNP positions (CHROM, POS), one row per genotype column
#   variables.csv      one row per sample, same order as the genotype matrix, with the nine
#                      GWAS variables already transformed (written by analysis/04)
# Output
#   out.rds            list(p, q, z, gif, loci): genomic-control P-values, Q-values and
#                      Z-scores (SNPs x variables), inflation factors, SNP positions

suppressPackageStartupMessages({
  library(LEA)
  library(data.table)
  library(parallel)
  library(qvalue)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
  stop("Usage: Rscript gwas_lfmm.R <genotypes.lfmm> <positions.012.pos> <variables.csv> <out.rds> [chunk_size] [K] [n_cores]")
f_geno     <- args[1]
f_pos      <- args[2]
f_vars     <- args[3]
f_out      <- args[4]
chunk_size <- if (length(args) >= 5) as.integer(args[5]) else 200000L
K          <- if (length(args) >= 6) as.integer(args[6]) else 1L
n_cores    <- if (length(args) >= 7) as.integer(args[7]) else 1L

GWAS_VARS <- c("VIZ_COUNT", "MATURE_ACORNS", "IMMAT_ACORNS", "ENLARGED_CUPS", "FLOWERS",
               "CANOPY_CLOSURE", "SPRING_PHENO", "SOIL_PC1", "MIDNOV_LAI")

# ── Explanatory matrix ───────────────────────────────────────────────────────
X <- as.matrix(read.csv(f_vars)[, GWAS_VARS])

# ── Genotypes (any residual missing calls, coded 9, are mean-imputed) ────────
Y    <- as.matrix(fread(f_geno, header = FALSE))
loci <- fread(f_pos, header = FALSE, col.names = c("CHROM", "POS"))
stopifnot(nrow(Y) == nrow(X), ncol(Y) == nrow(loci))
miss <- Y == 9
if (any(miss)) {
  Y[miss] <- NA
  mu <- colMeans(Y, na.rm = TRUE)
  Y[miss] <- round(mu[col(Y)[miss]])
}
cat("Genotypes:", nrow(Y), "samples x", ncol(Y), "SNPs;", sum(miss), "missing calls imputed\n")

# ── LFMM2 fit (ridge estimator with K latent factors) ────────────────────────
fit <- lfmm2(input = Y, env = X, K = K)

# ── Association tests in SNP chunks; genomic control applied genome-wide ─────
starts <- seq(1L, ncol(Y), by = chunk_size)
ends   <- pmin(starts + chunk_size - 1L, ncol(Y))
test_chunk <- function(i)
  lfmm2.test(fit, input = Y[, starts[i]:ends[i], drop = FALSE], env = X,
             genomic.control = FALSE)$zscores
z_list <- if (n_cores > 1 && .Platform$OS.type == "unix")
  mclapply(seq_along(starts), test_chunk, mc.cores = n_cores) else
  lapply(seq_along(starts), test_chunk)
z <- do.call(cbind, z_list)                     # variables x SNPs
rm(Y, z_list); invisible(gc())

chi2 <- z^2
gif  <- apply(chi2, 1, median, na.rm = TRUE) / qchisq(0.5, df = 1, lower.tail = FALSE)
p    <- pchisq(chi2 / gif, df = 1, lower.tail = FALSE)
q    <- t(apply(p, 1, function(x) qvalue(x)$qvalues))

print(data.frame(Trait = GWAS_VARS, GIF = round(gif, 4)))
print(setNames(rowSums(q < 0.05), GWAS_VARS))

# Only Z-scores are stored (xz-compressed); P- and Q-values are recomputed
# exactly from them (analysis/04_genomic_associations.Rmd, lfmm_pq()).
z <- t(z); colnames(z) <- GWAS_VARS
saveRDS(list(z = z, gif = data.frame(Trait = GWAS_VARS, GIF = round(gif, 4)), loci = loci),
        f_out, compress = "xz")
cat("Saved:", f_out, "\n")
