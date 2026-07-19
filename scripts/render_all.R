library(rmarkdown)

# Set pandoc path if not on system PATH (RStudio installation)
if (nchar(Sys.getenv("RSTUDIO_PANDOC")) == 0) {
  candidates <- c(
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/quarto/bin/tools"
  )
  for (p in candidates) {
    if (dir.exists(p)) { Sys.setenv(RSTUDIO_PANDOC = p); break }
  }
}

rmd_dir <- "analysis"
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rmds <- c(
  "01_descriptive_analysis.Rmd",
  "02_microclimate_analysis.Rmd",
  "03_statistical_models.Rmd"
)

for (f in rmds) {
  rmd_path <- file.path(rmd_dir, f)
  cat("\n══ Rendering:", f, "══\n")
  tryCatch(
    rmarkdown::render(
      input      = rmd_path,
      output_dir = out_dir,
      envir      = new.env()
    ),
    error = function(e) cat("✗ ERROR in", f, ":\n", conditionMessage(e), "\n")
  )
  cat("✓ Done:", f, "\n")
}

cat("\nAll done. HTML files in:", out_dir, "\n")
