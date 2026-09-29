# Environmental and genomic associations of inter-annual and inter-tree variation in acorn production in pedunculate oak

Data and code for the analyses of acorn production by 38 pedunculate oak
(*Quercus robur*) trees at Wytham Woods, Oxford, UK, monitored from 2020 to
2023.

## Overview

Oaks mast: acorn crops vary strongly from year to year and are synchronised
across a population, yet individual trees within the same woodland differ
widely in how many acorns they produce. This project asks why some trees are
consistent "super-producers" while others contribute little, and at which
stage of acorn development the difference arises.

Mature acorns were counted on 38 trees each year from 2020 to 2023, and
litter traps under each crown collected the earlier reproductive structures
that trees shed (flowers, enlarged cups, immature acorns) together with
knopper galls, from 2021 to 2023. These counts are related to tree size,
canopy closure, soil chemistry and texture, canopy phenology, trunk
microclimate and powdery mildew load, and to genome-wide SNP variation from
whole-genome sequencing of 37 of the trees.

The analyses:

1. describe inter-annual and inter-tree variation in acorn production,
   identify high-, medium- and low-producing trees and test how robust these
   classes are;
2. summarise weather and trunk microclimate over the study period and test
   whether microclimate differs among producer classes;
3. select environmental and tree-level predictors of each reproductive stage
   with random forests (VSURF), interpret them with Shapley values, and test
   them with Gamma generalised linear models and mixed models;
4. run an exploratory genome-wide association analysis (LFMM2) of the
   reproductive traits and their significant environmental predictors, with
   SNP annotation and gene-set enrichment.

## Workflow

The workflow runs end to end from the raw data in `data/raw/` to every table
and figure in the article. Four R Markdown notebooks follow the four
subsections of the Results. Each notebook is rendered to GitHub Markdown and
embeds the figures it produces.

| Notebook | Results subsection | Main text | Supplementary Information |
|---|---|---|---|
| [`01_acorn_production`](analysis/01_acorn_production.md) | Inter-annual variability in acorn production | Figures 1–3 | Tables 1–5, Figure 11 |
| [`02_weather_microclimate`](analysis/02_weather_microclimate.md) | Inter-annual variability in weather | – | Figures 1–4, Table 10 |
| [`03_environmental_associations`](analysis/03_environmental_associations.md) | Environmental and tree-level associates of acorn production | Tables 1–2, Figure 4 | Figures 5–10, 12; Tables 6–9, 11 |
| [`04_genomic_associations`](analysis/04_genomic_associations.md) | Genomic associates of acorn production | Figure 5 | Figures 13–17; Tables 12–16 |

Notebooks 2–4 use the tree-level tables that notebook 1 writes to
`data/derived/`, so they must be run in order.

## Repository structure

```
.
├── analysis/            R Markdown notebooks (.Rmd) and rendered Markdown (.md)
├── R/setup.R            shared paths, study-tree panel, plot theme, data loaders
├── data/
│   ├── raw/
│   │   ├── reproduction/   visual acorn counts, litter-trap counts, powdery mildew
│   │   ├── trees/          tree, soil and canopy-phenology attributes
│   │   ├── climate/        tree-mounted logger temperatures, ECN rainfall
│   │   └── genomics/       genotypes, PCA, kinship, sequencing statistics, reference annotation
│   └── derived/         tree-level tables written by the notebooks
├── cache/               stored results of computationally intensive steps
│   ├── vsurf/           VSURF random-forest variable selection
│   ├── shapley/         Shapley values
│   └── gwas/            sNMF, LFMM2 association results, SNP annotation, gene-set enrichment
├── genomics/            scripts run on a computing cluster (sNMF, imputation, LFMM2)
├── outputs/
│   ├── figures/         main and supplementary figures (PNG and PDF); multi-panel
│   │                    figures are saved as one file per panel
│   └── tables/          main and supplementary tables (CSV) and key numbers per notebook
└── render_all.R         renders the four notebooks in order
```

## Requirements

The analyses were run with R 4.5.2 on Windows; any recent R (≥ 4.3) should
work. Rendering the notebooks needs pandoc, which is bundled with RStudio.

| Source | Packages |
|---|---|
| CRAN | `tidyverse`, `readxl`, `patchwork`, `lubridate`, `ggrepel`, `ggsci`, `ggcorrplot`, `VSURF`, `ranger`, `iml`, `car`, `lme4`, `broom.mixed`, `multcomp`, `emmeans`, `DHARMa`, `MASS`, `vegan`, `data.table`, `rmarkdown` |
| Bioconductor | `LEA`, `qvalue`, `GenomicRanges`, `rtracklayer`, `clusterProfiler`, `org.At.tair.db` |

```r
install.packages(c("tidyverse", "readxl", "patchwork", "lubridate", "ggrepel",
                   "ggsci", "ggcorrplot", "VSURF", "ranger", "iml", "car",
                   "lme4", "broom.mixed", "multcomp", "emmeans", "DHARMa",
                   "MASS", "vegan", "data.table", "rmarkdown", "BiocManager"))
BiocManager::install(c("LEA", "qvalue", "GenomicRanges", "rtracklayer",
                       "clusterProfiler", "org.At.tair.db"))
```

Gene-set enrichment of KEGG pathways downloads the current KEGG annotation,
so it needs an internet connection the first time it is computed. The exact
package versions used are listed in the `sessionInfo()` at the end of each
rendered notebook.

## Reproducing the analyses

Open `wytham-oak-masting.Rproj` (or set the working directory to the
repository root) and run

```r
source("render_all.R")          # all notebooks
```

or `Rscript render_all.R 03` for a single notebook. A full run from the
cached results takes a few minutes.

Steps whose results are stored in `cache/` (Shapley values, LFMM2, SNP
annotation, gene-set enrichment) are read from the cache; deleting a cache
file makes the notebook recompute it. The VSURF variable selections are read
from the fitted objects in `cache/vsurf/`: random-forest variable selection
depends on the random number stream, which differs between package versions,
so a re-run may retain a slightly different set of variables. Read
processing, variant calling, filtering, LD pruning, PCA, kinship estimation,
sNMF and genotype imputation were run on a computing cluster and are
described in [`genomics/README.md`](genomics/README.md).

## Data

| File (in `data/raw/`) | Content |
|---|---|
| `reproduction/visual_counts_2020_2023.xlsx` | Mature acorns counted by visual survey, per tree and year |
| `reproduction/litter_traps_2021_2023.xlsx` | Litter-trap contents (mature and immature acorns, enlarged cups, flowers, galls), mean per 0.25 m² trap, per tree and year |
| `reproduction/powdery_mildew_2021.xlsx` | Dry weight of abscised, mildew-infected leaves per trap, August 2021 |
| `trees/tree_attributes.xlsx` | Location, tree dimensions, canopy, soil chemistry and texture, canopy phenology (LAI), summed kinship |
| `climate/tree_microclimate_daily_2020_2023.rds` | Daily temperature from trunk-mounted loggers |
| `climate/wytham_daily_precipitation_2020_2023.csv` | Daily rainfall, UK Environmental Change Network (ECN) Wytham station |
| `genomics/` | See [`genomics/README.md`](genomics/README.md) |

The raw files also contain records for two trees outside the 38-tree study
panel (Tree 10, with no records after 2021, and an extra tree, ET, counted in
2023 only). These are removed when the data are loaded (`STUDY_TREES` in
`R/setup.R`).

Raw sequencing reads are available from the NCBI Sequence Read Archive (see
the article's data availability statement).

## Citation

McClory R., Ellis R. H., Hung T. H., Lukac M., Clark J., Sheldon B. C.,
MacKay J. J. Environmental and genomic associations of inter-annual and
inter-tree variation in acorn production in pedunculate oak.
*Forestry: An International Journal of Forest Research*.

## Contact

Questions about the code and data: a@hung.bio
