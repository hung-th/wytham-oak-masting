# Render all analysis notebooks to GitHub-flavoured Markdown, in order.
# Run from the project root:  Rscript render_all.R        (all notebooks)
#                             Rscript render_all.R 03     (a single notebook)

if (!nzchar(Sys.getenv("RSTUDIO_PANDOC")) && !rmarkdown::pandoc_available()) {
  for (p in c("C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
              "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools",
              "/usr/lib/rstudio/resources/app/bin/quarto/bin/tools"))
    if (dir.exists(p)) { Sys.setenv(RSTUDIO_PANDOC = p); break }
}

notebooks <- c(
  "01_acorn_production.Rmd",
  "02_weather_microclimate.Rmd",
  "03_environmental_associations.Rmd",
  "04_genomic_associations.Rmd"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0)
  notebooks <- notebooks[substr(notebooks, 1, 2) %in% args]

for (nb in notebooks) {
  message("\n== Rendering ", nb, " ==")
  rmarkdown::render(file.path("analysis", nb), envir = new.env(), quiet = TRUE)
}
message("\nDone.")
