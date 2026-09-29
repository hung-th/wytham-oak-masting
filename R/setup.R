# Shared configuration for all analysis notebooks.
# Sourced at the top of every analysis/*.Rmd (knitted with the project root as
# the working directory, see render_all.R).

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(patchwork)
})

# ── Paths ────────────────────────────────────────────────────────────────────
DIR_RAW     <- file.path("data", "raw")
DIR_DERIVED <- file.path("data", "derived")
DIR_CACHE   <- "cache"
DIR_FIG     <- file.path("outputs", "figures")
DIR_TAB     <- file.path("outputs", "tables")
for (d in c(DIR_DERIVED, DIR_FIG, DIR_TAB,
            file.path(DIR_CACHE, c("vsurf", "shapley", "gwas"))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Study panel ──────────────────────────────────────────────────────────────
# The study panel comprises 38 trees. The raw files also contain records for
# two further trees that are not part of the panel and are removed on load:
#   Tree 10 – no records after 2021
#   ET      – extra tree, visual count in 2023 only
STUDY_TREES <- as.character(c(1:6, 8, 9, 11:39, 41))
stopifnot(length(STUDY_TREES) == 38)

STUDY_YEARS <- 2020:2023

# ── Plot theme and palettes ──────────────────────────────────────────────────
proj_theme <- theme_classic(base_size = 11) +
  theme(
    axis.title       = element_text(size = 12, colour = "black"),
    axis.text        = element_text(size = 10, colour = "black"),
    axis.line        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.6, colour = "black"),
    legend.title     = element_text(size = 10, colour = "black"),
    legend.text      = element_text(size = 9,  colour = "black"),
    strip.text       = element_text(size = 11, colour = "black"),
    strip.background = element_blank(),
    strip.clip       = "off",
    panel.grid       = element_blank()
  )
theme_set(proj_theme)

YEAR_COLS  <- c("2020" = "#BC3C29", "2021" = "#7876B1",
                "2022" = "#0072B5", "2023" = "#20854E")
CLASS_COLS <- c("High" = "#BC3C29", "Medium" = "#0072B5", "Low" = "grey45")

# ── Output helpers ───────────────────────────────────────────────────────────
# Save a ggplot/patchwork object as PDF + PNG and embed the PNG in the
# rendered Markdown (image links are relative to analysis/*.md).
save_fig <- function(plot, name, width, height, dpi = 300, show = TRUE) {
  ggsave(file.path(DIR_FIG, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIG, paste0(name, ".png")), plot,
         width = width, height = height, units = "mm", dpi = dpi)
  if (show) show_fig(name)
  invisible(plot)
}

show_fig <- function(name) {
  cat(sprintf("\n![](../outputs/figures/%s.png)\n\n", name))
}

save_tab <- function(df, name) {
  write_csv(df, file.path(DIR_TAB, paste0(name, ".csv")), na = "")
  invisible(df)
}

fmt_p <- function(p, digits = 4) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = digits)))
}

sig_stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

# ── Raw-data loaders (study panel only) ──────────────────────────────────────
# Visual survey of mature acorns: one count per tree per year (2020–2023).
load_visual_counts <- function() {
  read_excel(file.path(DIR_RAW, "reproduction", "visual_counts_2020_2023.xlsx")) |>
    transmute(Tree  = as.character(Tree),
              Year  = as.integer(format(as.Date(Date), "%Y")),
              Count = as.numeric(Counts)) |>
    filter(Tree %in% STUDY_TREES)
}

# Litter traps: mean count per 0.25 m² trap per tree per year (2021–2023).
# All traps of Trees 3, 6, 32, 34 and 36 were destroyed in 2023 (no record).
load_litter_traps <- function() {
  read_excel(file.path(DIR_RAW, "reproduction", "litter_traps_2021_2023.xlsx")) |>
    transmute(Tree          = as.character(Tree),
              Year          = as.integer(Year),
              MATURE_ACORNS = Acorns,
              IMMAT_ACORNS  = Immature_Acorns,
              ENLARGED_CUPS = Enlarged_cup,
              FLOWERS       = Flower,
              GALLS         = Gall) |>
    filter(Tree %in% STUDY_TREES)
}

LITTER_VARS <- c("MATURE_ACORNS", "IMMAT_ACORNS", "ENLARGED_CUPS",
                 "FLOWERS", "GALLS")

# Readable labels used in figures and tables
VAR_LABELS <- c(
  VIZ_COUNT      = "Mature acorns (visual count)",
  MATURE_ACORNS  = "Mature acorns (litter trap)",
  IMMAT_ACORNS   = "Immature acorns",
  ENLARGED_CUPS  = "Enlarged cups",
  FLOWERS        = "Flowers",
  GALLS          = "Galls",
  ELEVATION      = "Site Elevation",
  ASPECT         = "Slope Aspect",
  HEIGHT         = "Tree Height",
  DBH            = "DBH",
  CROWN_SPREAD   = "Crown Spread",
  CROWN_CROSS    = "Crown Cross",
  CANOPY_CLOSURE = "Canopy Closure",
  STEM_DENSITY   = "Stem Density",
  PH             = "pH",
  P_MGL          = "Phosphorus",
  K_MGL          = "Potassium",
  Mg_MGL         = "Magnesium",
  SAND           = "Sand",
  SILT           = "Silt",
  CLAY           = "Clay",
  SOIL_CLASS     = "Soil Class",
  SPRING_PHENO   = "Spring Phenology",
  MIDAPRIL_LAI   = "Mid-April LAI",
  AUTUMN_PHENO   = "Autumn Phenology",
  MIDNOV_LAI     = "Mid-November LAI",
  MIDDEC_LAI     = "Mid-December LAI",
  AUPPC          = "AUPPC",
  SOIL_PC1       = "Soil PC1"
)
