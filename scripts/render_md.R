library(rmarkdown)

if (nchar(Sys.getenv("RSTUDIO_PANDOC")) == 0) {
  candidates <- c(
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/quarto/bin/tools"
  )
  for (p in candidates) {
    if (dir.exists(p)) { Sys.setenv(RSTUDIO_PANDOC = p); break }
  }
}

# ── Analysis Rmds → analysis/ (figures land in analysis/*_files/) ────────────
analysis_rmds <- c(
  "01_descriptive_analysis.Rmd",
  "02_microclimate_analysis.Rmd",
  "03_statistical_models.Rmd"
)

for (f in analysis_rmds) {
  rmd_path <- file.path("analysis", f)
  cat("\n== Rendering (MD):", f, "==\n")
  tryCatch(
    rmarkdown::render(
      input         = rmd_path,
      output_format = github_document(html_preview = FALSE),
      output_dir    = "analysis",
      envir         = new.env()
    ),
    error = function(e) cat("ERROR in", f, ":\n", conditionMessage(e), "\n")
  )
  cat("Done:", f, "\n")
}

# ── GWAS Rmd → gwas/ ─────────────────────────────────────────────────────────
cat("\n== Rendering (MD): gwas_analysis.Rmd ==\n")
tryCatch(
  rmarkdown::render(
    input         = "gwas/gwas_analysis.Rmd",
    output_format = github_document(html_preview = FALSE),
    output_dir    = "gwas",
    envir         = new.env()
  ),
  error = function(e) cat("ERROR in gwas_analysis.Rmd:\n", conditionMessage(e), "\n")
)
cat("Done: gwas_analysis.Rmd\n")

# ── Fix absolute Windows paths in MD src= attributes ─────────────────────────
fix_paths <- function(md_file, strip_prefix) {
  txt <- readLines(md_file, warn = FALSE)
  txt <- gsub(strip_prefix, "", txt, fixed = TRUE)
  writeLines(txt, md_file)
}

repo_root <- normalizePath(".", winslash = "/")
for (f in analysis_rmds) {
  md <- file.path("analysis", sub("\\.Rmd$", ".md", f))
  if (file.exists(md))
    fix_paths(md, paste0(repo_root, "/analysis/"))
}
fix_paths("gwas/gwas_analysis.md", paste0(repo_root, "/gwas/"))

cat("\nAll done. Commit analysis/*.md, analysis/*_files/, gwas/*.md, gwas/*_files/\n")
