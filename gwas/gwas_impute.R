#!/usr/bin/env Rscript
# gwas_impute.R  —  parallel sNMF + site-wise imputation
# Usage: Rscript gwas_impute.R <run_dir> [n_cores] [chunk_size]
# Example: Rscript gwas_impute.R n36 16 500000

suppressPackageStartupMessages({
  library(LEA)
  library(data.table)
  library(parallel)
})

# ── Arguments ────────────────────────────────────────────────────────────────
# Usage: Rscript gwas_impute.R <run_dir> [lfmm_file] [n_cores] [chunk_size]
# lfmm_file is optional; if omitted, auto-detected from run_dir.
args       <- commandArgs(trailingOnly = TRUE)
run_dir    <- if (length(args) >= 1) args[1] else "n37"

# Optional explicit lfmm path (detected if arg contains ".lfmm")
explicit_lfmm <- NULL
arg_offset    <- 1L
if (length(args) >= 2 && grepl("\\.lfmm$", args[2])) {
  explicit_lfmm <- args[2]
  arg_offset    <- 2L
}
n_cores    <- if (length(args) > arg_offset)     as.integer(args[arg_offset + 1L]) else max(1L, detectCores() - 1L)
chunk_size <- if (length(args) > arg_offset + 1L) as.integer(args[arg_offset + 2L]) else 500000L

# ── File paths ────────────────────────────────────────────────────────────────
if (!is.null(explicit_lfmm)) {
  prefix      <- gsub("\\.lfmm$", "", basename(explicit_lfmm))
  full_lfmm   <- explicit_lfmm
  # When an explicit file is given it serves as both full and pruned input
  pruned_lfmm <- explicit_lfmm
} else {
  lfmm_files  <- list.files(run_dir, pattern = "\\.lfmm$")
  lfmm_base   <- lfmm_files[!grepl("pruned|imputed|snmf_sub|_no9", lfmm_files)]
  prefix      <- gsub("\\.lfmm$", "", lfmm_base[1])
  full_lfmm   <- file.path(run_dir, paste0(prefix, ".lfmm"))
  pruned_lfmm <- file.path(run_dir, paste0(prefix, "_pruned.lfmm"))
}
imputed_lfmm <- file.path(run_dir, paste0(prefix, "_imputed.lfmm"))
snmf_lfmm    <- file.path(run_dir, paste0(prefix, "_snmf_sub.lfmm"))

cat("=== gwas_impute.R ===\n")
cat("Run dir    :", run_dir,      "\n")
cat("Prefix     :", prefix,       "\n")
cat("Cores      :", n_cores,      "\n")
cat("Chunk size :", chunk_size,   "\n\n")

# ── Step 1: sNMF on pruned subsampled matrix ──────────────────────────────────
snmf_proj_path <- gsub("\\.lfmm$", ".snmfProject", snmf_lfmm)

if (file.exists(snmf_proj_path)) {
  cat("[1/3] Reloading existing sNMF project...\n")
  proj <- load.snmfProject(snmf_proj_path)
} else {
  cat("[1/3] Reading pruned lfmm for sNMF subsample...\n")
  pruned <- as.matrix(fread(pruned_lfmm, header = FALSE))
  cat("      Pruned matrix:", nrow(pruned), "samples x", ncol(pruned), "SNPs\n")

  set.seed(225)
  sub_idx <- sample(ncol(pruned), min(100000L, ncol(pruned)))
  fwrite(as.data.table(pruned[, sub_idx]), file = snmf_lfmm,
         sep = "\t", col.names = FALSE)
  rm(pruned); gc()

  cat("      Running sNMF K=1..6 (3 reps)...\n")
  proj <- snmf(snmf_lfmm, K = 1:6, entropy = TRUE,
               CPU = n_cores, repetitions = 3, project = "new")
}

# Print cross-entropy and select best run for K=1
mean_ce <- sapply(1:6, function(k)
  mean(sapply(1:3, function(r)
    tryCatch(cross.entropy(proj, K = k, run = r), error = function(e) NA_real_))))
cat("      Cross-entropy by K:", paste(round(mean_ce, 4), collapse = "  "), "\n")
cat("      Min CE at K =", which.min(mean_ce), "(using K=1 for imputation)\n\n")

best_run <- which.min(sapply(1:3, function(r)
  tryCatch(cross.entropy(proj, K = 1, run = r), error = function(e) Inf)))
Q_mat <- Q(proj, K = 1, run = best_run)  # n_samples × 1
cat("      Q matrix extracted (", nrow(Q_mat), "samples x", ncol(Q_mat), "components)\n\n")

# ── Step 2: Parallel imputation ───────────────────────────────────────────────
if (file.exists(imputed_lfmm)) {
  cat("[2/3] Imputed file already exists — skipping.\n")
  cat("      Delete", imputed_lfmm, "to re-run.\n\n")
} else {
  cat("[2/3] Reading full genotype matrix...\n")
  Y <- as.matrix(fread(full_lfmm, header = FALSE))
  n_samples <- nrow(Y)
  n_snps    <- ncol(Y)
  n_missing <- sum(Y == 9L)
  cat("      Full matrix:", n_samples, "samples x", n_snps, "SNPs\n")
  cat("      Missing genotypes (coded 9):", n_missing,
      sprintf("(%.3f%%)\n\n", 100 * n_missing / (n_samples * n_snps)))

  # Core imputation function: imputes one chunk and writes it to a temp file.
  # Returns the temp file path. Q_mat is available via fork (unix) or clusterExport.
  impute_chunk_to_file <- function(i, col_start, col_end, Y, tmp_dir) {
    chunk <- Y[, col_start:col_end, drop = FALSE]
    K     <- ncol(Q_mat)

    for (j in seq_len(ncol(chunk))) {
      geno    <- chunk[, j]
      missing <- geno == 9L
      if (!any(missing)) next
      obs <- !missing

      # Estimate allele frequencies F_j by OLS from observed genotypes
      Q_obs <- Q_mat[obs, , drop = FALSE]
      y_obs <- geno[obs]
      QQ    <- crossprod(Q_obs)
      Qy    <- crossprod(Q_obs, y_obs)
      F_j   <- tryCatch(
        solve(QQ + diag(1e-6, K), Qy),     # tiny ridge for stability
        error = function(e) matrix(mean(y_obs) / 2, K, 1)
      )
      F_j <- pmax(0, pmin(1, F_j))         # clamp to [0, 1]

      # Sample Binomial(2, Q·F) for each missing individual ("random" method)
      p_miss <- Q_mat[missing, , drop = FALSE] %*% F_j
      p_miss <- pmax(0, pmin(1, p_miss))
      chunk[missing, j] <- rbinom(sum(missing), 2L, as.numeric(p_miss))
    }

    # Write this chunk (rows = samples, cols = SNP slice) to a temp file
    tmp_path <- file.path(tmp_dir, sprintf("chunk_%05d.tmp", i))
    data.table::fwrite(as.data.table(chunk), file = tmp_path,
                       sep = " ", col.names = FALSE)
    tmp_path
  }

  # Split SNP columns into chunks
  tmp_dir      <- file.path(run_dir, "impute_tmp")
  dir.create(tmp_dir, showWarnings = FALSE)
  chunk_starts <- seq(1L, n_snps, by = chunk_size)
  chunk_ends   <- pmin(chunk_starts + chunk_size - 1L, n_snps)
  n_chunks     <- length(chunk_starts)
  cat("      Splitting into", n_chunks, "chunks of up to", chunk_size, "SNPs\n")
  cat("      Temp files in:", tmp_dir, "\n")
  cat("      Using", n_cores, "parallel cores\n\n")

  # Run in parallel — each worker writes its chunk to its own temp file
  if (.Platform$OS.type == "unix") {
    cat("      [unix] Using mclapply (fork)\n")
    tmp_files <- mclapply(
      seq_len(n_chunks),
      function(i) {
        if (i %% 10 == 0)
          message(sprintf("        chunk %d / %d", i, n_chunks))
        impute_chunk_to_file(i, chunk_starts[i], chunk_ends[i], Y, tmp_dir)
      },
      mc.cores = n_cores
    )
  } else {
    cat("      [windows] Using socket cluster\n")
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("Q_mat", "impute_chunk_to_file", "tmp_dir",
                        "chunk_starts", "chunk_ends", "Y"), envir = environment())
    clusterEvalQ(cl, { library(data.table); set.seed(225) })
    tmp_files <- parLapply(
      cl,
      seq_len(n_chunks),
      function(i) impute_chunk_to_file(i, chunk_starts[i], chunk_ends[i], Y, tmp_dir)
    )
    stopCluster(cl)
  }
  tmp_files <- unlist(tmp_files)

  # Check all chunks completed
  missing_chunks <- tmp_files[!file.exists(tmp_files)]
  if (length(missing_chunks) > 0)
    stop("Some chunks failed: ", paste(missing_chunks, collapse = ", "))

  # ── Stitch temp files horizontally ─────────────────────────────────────────
  # R has a hard limit of 128 open connections. Read chunk files in batches of
  # 60, loading all n_samples lines per chunk into memory, then write row by row.
  # Peak RAM ≈ total size of all chunk files as strings (~600 MB for 8.5M SNPs).
  cat("\n[3/3] Stitching", n_chunks, "chunk files into", imputed_lfmm, "...\n")
  # R's hard connection limit is 128. Query how many are already open (stdin/stdout/
  # stderr + R internals), leave a buffer of 5, and cap at that.
  batch_size <- min(n_chunks, max(1L, 128L - nrow(showConnections()) - 5L))
  cat("      Connection batch size:", batch_size,
      "(128 limit -", nrow(showConnections()), "open -5 buffer)\n")
  all_lines  <- vector("list", n_chunks)

  for (b in seq(1L, n_chunks, by = batch_size)) {
    batch_idx <- b:min(b + batch_size - 1L, n_chunks)
    cat(sprintf("      reading chunks %d-%d / %d\n", b, max(batch_idx), n_chunks))
    in_cons <- lapply(tmp_files[batch_idx], file, open = "r")
    for (k in seq_along(batch_idx)) {
      all_lines[[batch_idx[k]]] <- readLines(in_cons[[k]])
    }
    lapply(in_cons, close)
  }

  out_con <- file(imputed_lfmm, open = "w")
  for (row_i in seq_len(n_samples)) {
    writeLines(paste(vapply(all_lines, `[`, character(1L), row_i), collapse = " "),
               out_con)
  }
  close(out_con)
  rm(all_lines); gc()

  # Clean up temp files
  unlink(tmp_files)
  unlink(tmp_dir, recursive = TRUE)
  cat("      Done. Temp files removed.\n\n")
}

cat("=== Complete ===\n")
cat("Download these files to continue the analysis locally:\n")
cat(" ", imputed_lfmm, "\n")
cat(" ", snmf_lfmm,    "\n")
cat(" ", snmf_proj_path, " + all K*/ subdirectories\n")
cat("  n36/pca.eigenvec\n")
cat("  n36/pca.eigenval\n")
cat("  n36/relatedness.relatedness2\n")
cat("  n36/", prefix, ".012.pos\n", sep = "")
