# Wytham Woods Oak Masting

Reproducible analysis code for the study of inter-individual variation in acorn production and its environmental and morphological predictors in *Quercus robur* at Wytham Woods, Oxford.

## Repository structure

```
wytham-oak-masting/
├── analysis/               R Markdown analysis documents
│   ├── 01_descriptive_analysis.Rmd     Descriptive statistics and annual production figures
│   ├── 02_microclimate_analysis.Rmd    Microclimate characterisation and visualisation
│   ├── 03_statistical_models.Rmd       VSURF variable selection, GLMs, GLMMs, Shapley values
│   └── proj_theme.R                    Shared ggplot2 theme
├── gwas/                   GWAS pipeline (LFMM2, sNMF, enrichment)
│   ├── gwas_analysis.Rmd               End-to-end GWAS report
│   ├── gwas_impute.R                   Genotype imputation
│   ├── gwas_lfmm.R                     LFMM2 association tests
│   ├── run_snmf.R                      Population structure (sNMF)
│   └── main.R                          Orchestration script
├── data/
│   ├── raw/                Input data (Excel / CSV / RDS)
│   └── processed/          VSURF and Shapley cached results (.rds)
├── outputs/
│   ├── figures/            Publication figures (PDF + PNG)
│   └── tables/             Result tables (CSV)
└── scripts/
    └── render_all.R        Render all three analysis Rmds to HTML
```

## Reproducing the analysis

```r
# From the repo root
source("scripts/render_all.R")
```

This renders `01_`, `02_`, and `03_` in order. VSURF and Shapley computations are computationally intensive; cached `.rds` files in `data/processed/` are used by default (`eval_vsurf <- FALSE` in 03). Set `eval_vsurf <- TRUE` to re-run from scratch.

## Data

Raw data are in `data/raw/`. Key files:

| File | Contents |
|---|---|
| `wytham_full_dataset.xlsx` | Tree-level predictors and mean annual responses |
| `litter_traps_2021_2023.xlsx` | Annual litter trap counts per tree |
| `visual_counts_2020_2023.xlsx` | Annual visual acorn counts |
| `tree_microclimate_daily_2020_2023.rds` | Daily temperature logger data |
| `wytham_daily_precipitation_2020_2023.csv` | External weather station rainfall |

## GWAS data

Large genomic files (`.lfmm`, `.012`, `.snmf/`) are excluded from git due to size (see `.gitignore`). They should be obtained from the associated data archive (Zenodo/OSF). Place them in `gwas/` alongside the `.Rmd` before rendering `gwas_analysis.Rmd`.

## Dependencies

All analyses run in R ≥ 4.3. Key packages: `tidyverse`, `readxl`, `VSURF`, `iml`, `lme4`, `lmerTest`, `broom.mixed`, `car`, `DHARMa`, `patchwork`, `ggcorrplot`. GWAS additionally requires `LEA`, `data.table`, `qvalue`, `vegan`.
