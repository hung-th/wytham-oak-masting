Wytham Masting Study — Descriptive Analysis
================
Wytham Masting Project
2026-07-29

# Overview

This document reproduces all **descriptive analyses** for the manuscript
*“Environmental and genomic associations of inter-annual and inter-tree
variation in acorn production in pedunculate oak”*.

Figures produced here correspond to:

| Figure    | Description                                              |
|-----------|----------------------------------------------------------|
| Fig. 1a   | Annual acorn counts — visual method, per tree            |
| Fig. 1b–f | Annual reproductive material — litter traps              |
| Fig. 2    | % contribution per tree to total annual production       |
| Fig. 3    | Correlation: visual counts vs. litter trap mature acorns |
| Table S1  | Descriptive statistics for all response variables        |
| Table S2  | Population-level variability, synchrony, autocorrelation |

------------------------------------------------------------------------

# Packages and paths

``` r
library(tidyverse)     # data manipulation and ggplot2
library(readxl)        # read Excel files
library(ggthemes)      # extra ggplot2 themes
library(viridis)       # colour palettes
library(geomtextpath)  # smoothed labelled lines (Fig. 3)
library(plotrix)       # std.error()
library(car)           # Anova() for type-II tests
library(scales)        # scale formatting
library(patchwork)     # multi-panel figures

# ── Paths (relative to R1_revision/) ──────────────────────────────────────
DATA_RAW <- file.path("..", "data", "raw")
FIG_OUT  <- file.path("..", "outputs", "figures")
TAB_OUT  <- file.path("..", "outputs", "tables")

dir.create(FIG_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_OUT, recursive = TRUE, showWarnings = FALSE)

source("proj_theme.R")
select <- dplyr::select
filter <- dplyr::filter
```

------------------------------------------------------------------------

# Data loading

## Visual counts (mature acorn counts, survey method)

Visual counts were made annually by counting acorns on each tree via
binoculars. Data span 2020–2023. Tree “ET” (extra tree) is excluded from
main analyses as it was not part of the core 38-tree panel.

``` r
viz_raw <- read_excel(file.path(DATA_RAW, "visual_counts_2020_2023.xlsx"),
                      sheet = "Sheet1")

# Parse year from Date; remove ET
viz <- viz_raw %>%
  mutate(year = as.factor(format(as.Date(Date), "%Y")),
         Tree = as.factor(Tree)) %>%
  filter(Tree != "ET")

cat("Visual count data: ", nrow(viz), "rows |",
    "Years:", paste(levels(viz$year), collapse = ", "), "\n")
```

    ## Visual count data:  158 rows | Years: 2020, 2021, 2022, 2023

``` r
glimpse(viz)
```

    ## Rows: 158
    ## Columns: 4
    ## $ Date   <dttm> 2020-09-04, 2020-09-04, 2020-09-04, 2020-09-04, 2020-09-04, 20…
    ## $ Tree   <fct> 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, …
    ## $ Counts <dbl> 254, 76, 163, 146, 53, 39, 34, 92, 9, 0, 214, 118, 312, 162, 28…
    ## $ year   <fct> 2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020, 2020, 202…

## Litter trap data

Litter traps (0.25 m² per trap, six traps per tree) were deployed from
2021 onwards. **No litter trap data exist for 2020** — traps were not in
place during that season. This explains the absence of 2020 data in Fig.
1b–f.

``` r
trap <- read_excel(file.path(DATA_RAW, "litter_traps_2021_2023.xlsx"),
                   sheet = "Sheet1") %>%
  mutate(Year = as.factor(Year),
         Tree = as.factor(Tree))

cat("Litter trap data:", nrow(trap), "rows |",
    "Years:", paste(levels(trap$Year), collapse = ", "), "\n")
```

    ## Litter trap data: 115 rows | Years: 2021, 2022, 2023

``` r
glimpse(trap)
```

    ## Rows: 115
    ## Columns: 7
    ## $ Tree            <fct> 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17,…
    ## $ Year            <fct> 2021, 2021, 2021, 2021, 2021, 2021, 2021, 2021, 2021, …
    ## $ Acorns          <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ Immature_Acorns <dbl> 0.500000, 0.250000, 0.250000, 0.500000, 0.000000, 0.00…
    ## $ Enlarged_cup    <dbl> 0.7500000, 0.5000000, 1.0000000, 0.2500000, 1.0000000,…
    ## $ Flower          <dbl> 64.00000, 29.00000, 28.50000, 31.50000, 35.33333, 28.0…
    ## $ Gall            <dbl> 7.7500000, 3.7500000, 1.7500000, 2.5000000, 0.3333333,…

## Full dataset (tree-level explanatory + response variables)

This is the master dataset used for statistical modelling (Section 03).
Loaded here for descriptive statistics.

``` r
full_df <- read_excel(file.path(DATA_RAW, "wytham_full_dataset.xlsx"),
                      sheet = "Data") %>%
  filter(TREE_ID != "ET") %>%
  mutate(
    DBH           = as.numeric(DBH),
    VIZ_COUNT     = as.numeric(VIZ_COUNT),   # stored as character in Excel
    PRODUCER      = as.factor(PRODUCER),
    SOIL_CLASS    = as.factor(SOIL_CLASS),
    across(c(CLAY, SILT, SAND, Mg_MGL, K_MGL, P_MGL, PH), as.numeric)
  ) %>%
  # Remove nutrient index columns (retained as mg/L concentrations only)
  dplyr::select(-any_of(c("K_INDEX", "P_INDEX", "Mg_INDEX")))

cat("Full dataset:", nrow(full_df), "rows |",
    nlevels(as.factor(full_df$TREE_ID)), "trees\n")
```

    ## Full dataset: 39 rows | 39 trees

------------------------------------------------------------------------

# Descriptive statistics table for response variables

Response variables defined:

| Variable       | Definition                                     | Units |
|----------------|------------------------------------------------|-------|
| VIZ_COUNT      | Mature acorn count, visual method              | count |
| MATURE_ACORNS  | Mature acorns in litter traps (per 0.25 m²)    | count |
| IMMAT_ACORNS   | Immature acorns in litter traps                | count |
| ABORTED_ACORNS | Enlarged cups (aborted acorns) in litter traps | count |
| FLOWERS        | Catkin/flower counts in litter traps           | count |
| CROP_PERCENT   | % of total annual acorn crop                   | %     |

``` r
resp_vars <- c("VIZ_COUNT", "MATURE_ACORNS", "IMMAT_ACORNS",
               "ABORTED_ACORNS", "FLOWERS", "CROP_PERCENT")

desc_stats <- full_df %>%
  # Use one row per tree (the full_df may be tree-averaged already)
  distinct(TREE_ID, .keep_all = TRUE) %>%
  select(all_of(resp_vars)) %>%
  pivot_longer(everything(), names_to = "Variable") %>%
  group_by(Variable) %>%
  summarise(
    N      = sum(!is.na(value)),
    Mean   = mean(value, na.rm = TRUE),
    SD     = sd(value, na.rm = TRUE),
    Median = median(value, na.rm = TRUE),
    Min    = min(value, na.rm = TRUE),
    Max    = max(value, na.rm = TRUE),
    CV_pct = SD / Mean * 100,
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  arrange(match(Variable, resp_vars))

# Rename for readability
desc_stats_display <- desc_stats %>%
  rename(
    `Variable`   = Variable,
    `N`          = N,
    `Mean`       = Mean,
    `SD`         = SD,
    `Median`     = Median,
    `Min`        = Min,
    `Max`        = Max,
    `CV (%)`     = CV_pct
  )

knitr::kable(desc_stats_display,
             caption = "**Table S1.** Descriptive statistics for response
             variables measured in 38 pedunculate oak trees at Wytham Woods.
             CV = coefficient of variation (SD/Mean × 100). Note that
             MATURE_ACORNS, IMMAT_ACORNS, ABORTED_ACORNS and FLOWERS are
             averaged across litter trap collections within year. The high CVs
             justify the use of a Gamma (log link) GLM rather than Gaussian
             models (Gamma is appropriate for skewed, strictly positive
             continuous data).",
             align = "lrrrrrrr")
```

| Variable       |   N |  Mean |    SD | Median |  Min |    Max | CV (%) |
|:---------------|----:|------:|------:|-------:|-----:|-------:|-------:|
| VIZ_COUNT      |  39 | 65.03 | 64.84 |  38.67 | 1.67 | 238.00 |  99.71 |
| MATURE_ACORNS  |  39 |  2.96 |  4.76 |   1.11 | 0.00 |  24.89 | 160.96 |
| IMMAT_ACORNS   |  39 |  4.61 |  4.52 |   2.94 | 0.00 |  15.19 |  98.07 |
| ABORTED_ACORNS |  39 |  5.27 |  4.30 |   4.39 | 0.00 |  17.53 |  81.57 |
| FLOWERS        |  39 | 39.31 | 21.13 |  37.14 | 4.72 |  94.72 |  53.76 |
| CROP_PERCENT   |  39 |  0.03 |  0.03 |   0.02 | 0.00 |   0.10 |  98.51 |

**Table S1.** Descriptive statistics for response variables measured in
38 pedunculate oak trees at Wytham Woods. CV = coefficient of variation
(SD/Mean × 100). Note that MATURE_ACORNS, IMMAT_ACORNS, ABORTED_ACORNS
and FLOWERS are averaged across litter trap collections within year. The
high CVs justify the use of a Gamma (log link) GLM rather than Gaussian
models (Gamma is appropriate for skewed, strictly positive continuous
data).

``` r
# Save to outputs
write.csv(desc_stats_display,
          file.path(TAB_OUT, "Table_S1_descriptive_statistics_response_vars.csv"),
          row.names = FALSE)
```

------------------------------------------------------------------------

# Figure 1: Annual variation in acorn production

## Fig. 1a — Visual counts by year

Trees 27 (consistently poor producer) and 36 (consistently
super-producer) are marked with distinct point shapes and labelled.

``` r
# Tag focal trees
viz_tagged <- viz %>%
  mutate(focal = case_when(
    Tree == "36" ~ "Tree 36 (super producer)",
    Tree == "27" ~ "Tree 27 (poor producer)",
    TRUE         ~ "Other trees"
  ),
  focal = factor(focal, levels = c("Tree 36 (super producer)",
                                   "Tree 27 (poor producer)",
                                   "Other trees")))

fig1a <- ggplot(viz_tagged, aes(x = year, y = Counts)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", colour = "grey40",
               width = 0.5) +
  geom_jitter(data = filter(viz_tagged, focal == "Other trees"),
              colour = "grey60", size = 1.5, alpha = 0.7,
              width = 0.15, height = 0) +
  geom_point(data = filter(viz_tagged, focal != "Other trees"),
             aes(shape = focal, colour = focal),
             size = 2, stroke = 1.2) +
  scale_shape_manual(values = c(17, 15)) +
  scale_colour_manual(values = c("firebrick", "steelblue4")) +
  labs(x = "", y = "Acorn count (visual)",
       shape = "", colour = "") +
  guides(shape  = guide_legend(nrow = 2),
         colour = guide_legend(nrow = 2)) +
  proj_theme +
  theme(legend.position  = "top",
        legend.margin    = margin(t = -20),
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(10, "pt"),
        panel.grid = element_blank())

fig1a
```

<div class="figure">

<img src="C:/Users/th/Documents/wwo-masting/wytham-oak-masting/analysis/01_descriptive_analysis_files/figure-gfm/fig1a-1.png" alt="Fig. 1a. Acorn counts per tree per year at Wytham Woods (visual method). Grey boxplots show the population distribution; individual trees shown as jittered points. Trees 36 (super-producer) and 27 (poor producer) are highlighted." width="100%" />
<p class="caption">

Fig. 1a. Acorn counts per tree per year at Wytham Woods (visual method).
Grey boxplots show the population distribution; individual trees shown
as jittered points. Trees 36 (super-producer) and 27 (poor producer) are
highlighted.
</p>

</div>

``` r
ggsave(file.path(FIG_OUT, "Fig1a_visual_counts_by_year.pdf"),
       fig1a, width = 90, height = 100, units = "mm")
ggsave(file.path(FIG_OUT, "Fig1a_visual_counts_by_year.png"),
       fig1a, width = 60, height = 60, units = "mm", dpi = 300)
```

## Fig. 1b–f — Litter trap components by year

Litter traps were not deployed until 2021; the 2020 masting event is
captured only by the visual count method (Fig. 1a). This limitation is
explicitly stated in the figure caption.

``` r
# Reshape to long format for faceted plotting
trap_long <- trap %>%
  pivot_longer(cols = c(Acorns, Immature_Acorns, Enlarged_cup, Flower, Gall),
               names_to = "Category", values_to = "Count") %>%
  mutate(Category = factor(Category,
    levels = c("Acorns", "Immature_Acorns", "Enlarged_cup", "Flower", "Gall"),
    labels = c("Mature acorns", "Immature acorns",
               "Enlarged cups\n(aborted)", "Flowers",
               "Galls")))

fig1b_f <- ggplot(trap_long, aes(x = Year, y = Count)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", colour = "grey40",
               width = 0.5) +
  geom_jitter(colour = "black", size = 1.2, alpha = 0.7,
              width = 0.15, height = 0) +
  facet_wrap(~ Category, scales = "free_y", ncol = 3, axes = "all") +
  labs(x = "", y = "Count per trap (0.25 m²)") +
  proj_theme +
  theme(panel.grid = element_blank())

fig1b_f
```

<div class="figure">

<img src="C:/Users/th/Documents/wwo-masting/wytham-oak-masting/analysis/01_descriptive_analysis_files/figure-gfm/fig1b-f-panels-1.png" alt="Fig. 1b–f. Annual counts of reproductive and associated material from litter traps (0.25 m² per trap) under 38 pedunculate oak trees at Wytham Woods, 2021–2023. Note: litter traps were not deployed in 2020; that year's acorn production is captured by visual counts (Fig. 1a)." width="100%" />
<p class="caption">

Fig. 1b–f. Annual counts of reproductive and associated material from
litter traps (0.25 m² per trap) under 38 pedunculate oak trees at Wytham
Woods, 2021–2023. Note: litter traps were not deployed in 2020; that
year’s acorn production is captured by visual counts (Fig. 1a).
</p>

</div>

``` r
ggsave(file.path(FIG_OUT, "Fig1b-f_litter_trap_components.pdf"),
       fig1b_f, width = 170, height = 120, units = "mm")
ggsave(file.path(FIG_OUT, "Fig1b-f_litter_trap_components.png"),
       fig1b_f, width = 120, height = 120, units = "mm", dpi = 300)
```

------------------------------------------------------------------------

# Figure 2: Total acorn production per tree with annual deviations

Trees are ordered left to right by cumulative total count (highest
producer Tree 36 on left). Top panel: stacked raw counts by year. Lower
panels: deviation of each tree from the yearly population mean, one
panel per year (2021 excluded — no acorns matured).

``` r
# NEJM palette for years (2021 excluded from stacked bars — zero production)
nejm_4 <- c("2020" = "#BC3C29",
             "2022" = "#0072B5",
             "2023" = "#20854E")

# Order trees by cumulative total count
tree_order <- viz %>%
  group_by(Tree) %>%
  summarise(cum_count = sum(Counts, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(cum_count)) %>%
  pull(Tree)

viz_ordered <- viz %>%
  mutate(Tree = factor(Tree, levels = tree_order),
         year = as.character(year))

# ── Top panel: stacked raw counts (exclude 2021 — all zero) ─────────────
viz_nonzero <- viz_ordered %>% filter(year != "2021")

p_top <- ggplot(viz_nonzero, aes(x = Tree, y = Counts, fill = year)) +
  geom_col(width = 0.85, colour = NA) +
  scale_fill_manual(values = nejm_4, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Acorn counts") +
  guides(fill = guide_legend(nrow = 1, keywidth = 0.8, keyheight = 0.5)) +
  proj_theme +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = NA))

# ── Deviation panels: one per year ───────────────────────────────────────
deviation_panels <- function(yr) {
  d <- viz_ordered %>%
    filter(year == yr) %>%
    group_by(year) %>%
    mutate(deviation = Counts - mean(Counts, na.rm = TRUE)) %>%
    ungroup()

  ggplot(d, aes(x = Tree, y = deviation)) +
    geom_col(fill = "grey70", colour = "grey40", width = 0.8, linewidth = 0.2) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    annotate("text", x = Inf, y = Inf, label = yr,
             hjust = 1.1, vjust = 1.4, size = 3.5, fontface = "bold") +
    scale_y_continuous(limits = c(-150, 350)) +
    labs(x = NULL, y = NULL) +
    proj_theme +
    theme(axis.text.x = element_text(size = 7, colour = "black"))
}

p_2020 <- deviation_panels("2020") + theme(axis.text.x = element_blank(),
                                            axis.ticks.x = element_blank())
p_2022 <- deviation_panels("2022") +
  labs(y = "Deviation from yearly mean acorn count") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p_2023 <- deviation_panels("2023") +
  labs(x = "Tree")

fig2 <- p_top / p_2020 / p_2022 / p_2023 +
  plot_layout(heights = c(1.4, 1, 1, 1))

fig2
```

<div class="figure">

<img src="C:/Users/th/Documents/wwo-masting/wytham-oak-masting/analysis/01_descriptive_analysis_files/figure-gfm/fig2-counts-deviation-1.png" alt="Fig. 2. Total number of acorns produced by individual pedunculate oak trees at Wytham Woods over 2020, 2022, 2023 from visual counts (top); no acorns matured in 2021. Lower panels show deviation of each tree from the yearly mean acorn count. Trees are ranked from highest (Tree 36) to lowest (Tree 27) total production." width="100%" />
<p class="caption">

Fig. 2. Total number of acorns produced by individual pedunculate oak
trees at Wytham Woods over 2020, 2022, 2023 from visual counts (top); no
acorns matured in 2021. Lower panels show deviation of each tree from
the yearly mean acorn count. Trees are ranked from highest (Tree 36) to
lowest (Tree 27) total production.
</p>

</div>

``` r
ggsave(file.path(FIG_OUT, "Fig2_acorn_counts_and_deviations.pdf"),
       fig2, width = 150, height = 120, units = "mm")
ggsave(file.path(FIG_OUT, "Fig2_acorn_counts_and_deviations.png"),
       fig2, width = 150, height = 120, units = "mm", dpi = 300)

# Headline numbers
top7_trees <- tree_order[1:7]
top7_pct <- viz %>%
  mutate(year_total = ave(Counts, year, FUN = function(x) sum(x, na.rm=TRUE)),
         pct = Counts / year_total * 100) %>%
  filter(Tree %in% top7_trees) %>%
  group_by(year) %>%
  summarise(top7_pct = sum(pct, na.rm = TRUE))
cat("Top 7 trees % of annual total by year:\n")
```

    ## Top 7 trees % of annual total by year:

``` r
print(top7_pct)
```

    ## # A tibble: 4 × 2
    ##   year  top7_pct
    ##   <fct>    <dbl>
    ## 1 2020      38.8
    ## 2 2021       0  
    ## 3 2022      66.4
    ## 4 2023      58.3

------------------------------------------------------------------------

# Figure 3: Correlation between visual counts and litter trap acorns

Visual counts and litter trap mature acorn counts are compared to
validate the two methods. Trees 3, 6, 10, 32, 34, 36 are excluded as
they lack complete data across methods.

``` r
EXCLUDE_TREES <- c("3", "6", "10", "32", "34", "36")

viz_cor <- viz %>%
  filter(!Tree %in% EXCLUDE_TREES, year != "2020") %>%
  rename(Year = year)

trap_cor <- trap %>%
  filter(!Tree %in% EXCLUDE_TREES) %>%
  rename(Acorns_trap = Acorns)

combined <- inner_join(viz_cor, trap_cor, by = c("Tree", "Year"))

# Kendall and Pearson correlations reported in-text
kt <- cor.test(combined$Counts, combined$Acorns_trap, method = "kendall")
pr <- cor.test(log(combined$Counts + 1), log(combined$Acorns_trap + 1),
               method = "pearson")

cat("Kendall tau =", round(kt$estimate, 3), ", p =", round(kt$p.value, 4), "\n")
```

    ## Kendall tau = 0.718 , p = 0

``` r
cat("Pearson r (log+1) =", round(pr$estimate, 3), ", p =", round(pr$p.value, 4), "\n")
```

    ## Pearson r (log+1) = 0.783 , p = 0

``` r
nejm_2 <- c("2022" = "#0072B5", "2023" = "#20854E")

fig3 <- ggplot(combined, aes(x = Counts, y = Acorns_trap,
                              shape = Year, linetype = Year, colour = Year)) +
  geom_point(size = 2.5) +
  geom_smooth(aes(group = Year), method = "lm", formula = y ~ x,
              se = FALSE, linewidth = 0.8) +
  scale_shape_manual(values = c("2022" = 17, "2023" = 15)) +
  scale_linetype_manual(values = c("2022" = "solid",
                                   "2023" = "dashed")) +
  scale_colour_manual(values = nejm_2) +
  guides(shape    = guide_legend(title = NULL),
         linetype = guide_legend(title = NULL),
         colour   = guide_legend(title = NULL)) +
  labs(x = "Visual acorn count",
       y = "Mature acorns\nper trap (0.25 m²)") +
  proj_theme +
  theme(legend.position  = "top",
        legend.margin    = margin(t = 0, b = 0),
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(10, "pt"),
        legend.key.width  = unit(12, "pt"),
        panel.grid = element_blank())

fig3
```

<div class="figure">

<img src="C:/Users/th/Documents/wwo-masting/wytham-oak-masting/analysis/01_descriptive_analysis_files/figure-gfm/fig3-correlation-1.png" alt="Fig. 3. Correlation between visual acorn counts and litter trap mature acorn counts per 0.25 m² for individual trees, 2021–2023 (2020 excluded as litter traps not deployed). Lines fitted by ordinary least squares per year." width="100%" />
<p class="caption">

Fig. 3. Correlation between visual acorn counts and litter trap mature
acorn counts per 0.25 m² for individual trees, 2021–2023 (2020 excluded
as litter traps not deployed). Lines fitted by ordinary least squares
per year.
</p>

</div>

``` r
ggsave(file.path(FIG_OUT, "Fig3_visual_vs_trap_correlation.pdf"),
       fig3, width = 90, height = 90, units = "mm")
ggsave(file.path(FIG_OUT, "Fig3_visual_vs_trap_correlation.png"),
       fig3, width = 65, height = 60, units = "mm", dpi = 300)
```

------------------------------------------------------------------------

# Population-level reproductive variability, autocorrelation and synchrony

With only 4 years of data (3 non-zero years for litter traps), all
estimates of autocorrelation and synchrony carry considerable
uncertainty. Results should be interpreted as preliminary descriptors
rather than formal statistical tests.

``` r
# ── Individual-level CV ──────────────────────────────────────────────────
# Use NA for trees with mean_count = 0 (super-poor producers in all years)
indiv_cv <- viz %>%
  group_by(Tree) %>%
  summarise(
    mean_count = mean(Counts, na.rm = TRUE),
    sd_count   = sd(Counts, na.rm = TRUE),
    CV_pct     = ifelse(mean_count > 0, sd_count / mean_count * 100, NA_real_),
    n_years    = n(),
    .groups    = "drop"
  ) %>%
  arrange(desc(CV_pct))

cat("Individual-level CV (visual counts across all years):\n")
```

    ## Individual-level CV (visual counts across all years):

``` r
cat("  Mean CV:", round(mean(indiv_cv$CV_pct, na.rm = TRUE), 1), "%\n")
```

    ##   Mean CV: 127.3 %

``` r
cat("  Median CV:", round(median(indiv_cv$CV_pct, na.rm = TRUE), 1), "%\n")
```

    ##   Median CV: 117.8 %

``` r
cat("  Range:", round(min(indiv_cv$CV_pct, na.rm = TRUE), 1), "–",
    round(max(indiv_cv$CV_pct, na.rm = TRUE), 1), "%\n\n")
```

    ##   Range: 74.5 – 200 %

``` r
# ── Population CV per year ───────────────────────────────────────────────
# 2021 was a near-complete masting failure; if all counts = 0, CV is undefined
pop_cv <- viz %>%
  group_by(year) %>%
  summarise(
    mean_count = mean(Counts, na.rm = TRUE),
    sd_count   = sd(Counts, na.rm = TRUE),
    CV_pct     = ifelse(mean_count > 0, sd_count / mean_count * 100, NA_real_),
    .groups    = "drop"
  )

cat("Population CV by year (inter-individual variation):\n")
```

    ## Population CV by year (inter-individual variation):

``` r
print(pop_cv)
```

    ## # A tibble: 4 × 4
    ##   year  mean_count sd_count CV_pct
    ##   <fct>      <dbl>    <dbl>  <dbl>
    ## 1 2020       120.     107.    89.9
    ## 2 2021         0        0     NA  
    ## 3 2022        72.5    115.   159. 
    ## 4 2023        47.2     74.9  159.

``` r
cat("  Note: NA indicates a near-zero year where CV is undefined (mean ≈ 0).\n")
```

    ##   Note: NA indicates a near-zero year where CV is undefined (mean ≈ 0).

``` r
# ── Population-level variability (inter-annual CV of population mean) ────
pop_annual <- viz %>%
  group_by(year) %>%
  summarise(mean_count = mean(Counts), .groups = "drop")

cv_interannual <- sd(pop_annual$mean_count) / mean(pop_annual$mean_count) * 100
cat("\nInter-annual CV of population mean:", round(cv_interannual, 1), "%\n")
```

    ## 
    ## Inter-annual CV of population mean: 83.4 %

``` r
# ── Synchrony: Pearson correlation among all tree pairs ─────────────────
# Pivot to wide format (tree × year)
viz_wide <- viz %>%
  pivot_wider(id_cols = Tree, names_from = year, values_from = Counts)

# Compute pairwise Pearson correlations across years (requires ≥3 years)
sync_pairs <- combn(nrow(viz_wide), 2, function(idx) {
  x <- as.numeric(viz_wide[idx[1], -1])
  y <- as.numeric(viz_wide[idx[2], -1])
  # Use complete pairs only
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok], method = "pearson")
})

mean_sync <- mean(sync_pairs, na.rm = TRUE)
cat("\nPopulation synchrony (mean pairwise Pearson r across all tree pairs):\n")
```

    ## 
    ## Population synchrony (mean pairwise Pearson r across all tree pairs):

``` r
cat("  Mean r =", round(mean_sync, 3), "\n")
```

    ##   Mean r = 0.573

``` r
cat("  SD r =", round(sd(sync_pairs, na.rm = TRUE), 3), "\n")
```

    ##   SD r = 0.398

``` r
cat("  Note: based on only", sum(!is.na(sync_pairs)), "pairs with ≥3 shared years.\n")
```

    ##   Note: based on only 741 pairs with ≥3 shared years.

``` r
cat("  ⚠ With n = 4 years this estimate is unreliable; report with caution.\n\n")
```

    ##   ⚠ With n = 4 years this estimate is unreliable; report with caution.

``` r
# ── Lag-1 autocorrelation (per tree, then averaged) ──────────────────────
lag1_ac <- viz %>%
  arrange(Tree, year) %>%
  group_by(Tree) %>%
  summarise(
    ac1 = {
      x <- Counts
      # Require at least 3 time points and non-constant series
      if (length(x) >= 3 && sd(x, na.rm = TRUE) > 0 &&
          sd(x[-length(x)], na.rm = TRUE) > 0)
        suppressWarnings(cor(x[-length(x)], x[-1], method = "pearson"))
      else NA_real_
    },
    .groups = "drop"
  )

cat("Lag-1 autocorrelation of individual trees:\n")
```

    ## Lag-1 autocorrelation of individual trees:

``` r
cat("  Mean:", round(mean(lag1_ac$ac1, na.rm = TRUE), 3), "\n")
```

    ##   Mean: -0.718

``` r
cat("  SD:", round(sd(lag1_ac$ac1, na.rm = TRUE), 3), "\n")
```

    ##   SD: 0.232

``` r
cat("  ⚠ With n = 4 time points this is highly unreliable; treat descriptively.\n")
```

    ##   ⚠ With n = 4 time points this is highly unreliable; treat descriptively.

``` r
cat("  ⚠ The strongly negative mean value is dominated by the 2020 (high) → 2021\n")
```

    ##   ⚠ The strongly negative mean value is dominated by the 2020 (high) → 2021

``` r
cat("    (near-zero) transition. With only 3 lag-1 pairs per tree, this single\n")
```

    ##     (near-zero) transition. With only 3 lag-1 pairs per tree, this single

``` r
cat("    transition drives the estimate. Do not over-interpret the magnitude.\n")
```

    ##     transition drives the estimate. Do not over-interpret the magnitude.

``` r
# Summary table for manuscript (Table R2)
table_r2 <- tibble(
  Metric = c(
    "Mean individual CV (%)",
    "Population CV — 2020 (%)",
    "Population CV — 2021 (%) [mast failure]",
    "Population CV — 2022 (%)",
    "Population CV — 2023 (%)",
    "Inter-annual CV of population mean (%)",
    "Mean pairwise synchrony (Pearson r)",
    "Mean lag-1 autocorrelation (Pearson r)"
  ),
  Value = c(
    round(mean(indiv_cv$CV_pct, na.rm = TRUE), 1),
    round(filter(pop_cv, year == "2020")$CV_pct, 1),
    ifelse(is.na(filter(pop_cv, year == "2021")$CV_pct), "NA (mean ≈ 0)",
           as.character(round(filter(pop_cv, year == "2021")$CV_pct, 1))),
    round(filter(pop_cv, year == "2022")$CV_pct, 1),
    round(filter(pop_cv, year == "2023")$CV_pct, 1),
    round(cv_interannual, 1),
    round(mean_sync, 3),
    round(mean(lag1_ac$ac1, na.rm = TRUE), 3)
  ),
  Notes = c(
    "Mean across 38 trees; all years",
    "Inter-individual variation within year",
    "Inter-individual variation within year",
    "Inter-individual variation within year",
    "Inter-individual variation within year",
    "Year-to-year variation in population mean",
    "Mean of all pairwise r; n=4 yrs, use with caution",
    "Mean across trees; n=4 time points, very limited"
  )
)

knitr::kable(table_r2,
             caption = "**Table S2.** Reproductive variability and synchrony
             statistics for 38 pedunculate oak trees at Wytham Woods.
             CV = coefficient of variation. Synchrony estimated as mean
             pairwise Pearson correlation coefficient (Loreau & de Mazancourt
             2008 approach). Lag-1 autocorrelation reported per tree then
             averaged. All estimates based on n = 4 years of data and should
             be interpreted with caution.",
             align = "lrr")
```

| Metric | Value | Notes |
|:---|---:|---:|
| Mean individual CV (%) | 127.3 | Mean across 38 trees; all years |
| Population CV — 2020 (%) | 89.9 | Inter-individual variation within year |
| Population CV — 2021 (%) \[mast failure\] | NA (mean ≈ 0) | Inter-individual variation within year |
| Population CV — 2022 (%) | 158.9 | Inter-individual variation within year |
| Population CV — 2023 (%) | 158.8 | Inter-individual variation within year |
| Inter-annual CV of population mean (%) | 83.4 | Year-to-year variation in population mean |
| Mean pairwise synchrony (Pearson r) | 0.573 | Mean of all pairwise r; n=4 yrs, use with caution |
| Mean lag-1 autocorrelation (Pearson r) | -0.718 | Mean across trees; n=4 time points, very limited |

**Table S2.** Reproductive variability and synchrony statistics for 38
pedunculate oak trees at Wytham Woods. CV = coefficient of variation.
Synchrony estimated as mean pairwise Pearson correlation coefficient
(Loreau & de Mazancourt 2008 approach). Lag-1 autocorrelation reported
per tree then averaged. All estimates based on n = 4 years of data and
should be interpreted with caution.

``` r
write.csv(table_r2,
          file.path(TAB_OUT, "Table_S2_variability_synchrony.csv"),
          row.names = FALSE)
```

------------------------------------------------------------------------

# Super/medium/low producer classification excluding 2020

Producer categories are re-derived on 2021–2023 data to assess whether
the 2020 mast year drives the classification.

``` r
# Load full dataset PRODUCER labels (original classification on all 4 years)
full_prod <- read_excel(file.path(DATA_RAW, "wytham_full_dataset.xlsx"),
                        sheet = "Data") %>%
  filter(TREE_ID != "ET") %>%
  select(TREE_ID, PRODUCER) %>%
  distinct() %>%
  mutate(Tree = as.character(TREE_ID))

# Visual counts excluding 2020
viz_no2020 <- viz %>%
  filter(year != "2020") %>%
  group_by(Tree) %>%
  summarise(
    mean_count_no2020 = mean(Counts, na.rm = TRUE),
    sum_count_no2020  = sum(Counts,  na.rm = TRUE),
    .groups = "drop"
  )

# Re-classify using the same quantile thresholds applied to the no-2020 means
q33 <- quantile(viz_no2020$mean_count_no2020, 1/3, na.rm = TRUE)
q67 <- quantile(viz_no2020$mean_count_no2020, 2/3, na.rm = TRUE)

viz_no2020 <- viz_no2020 %>%
  mutate(
    Producer_no2020 = case_when(
      mean_count_no2020 >= q67 ~ "Super",
      mean_count_no2020 >= q33 ~ "Medium",
      TRUE                     ~ "Low"
    ),
    Producer_no2020 = factor(Producer_no2020,
                             levels = c("Super", "Medium", "Low"))
  )

# Merge with original classification
prod_compare <- full_prod %>%
  left_join(viz_no2020, by = "Tree") %>%
  mutate(PRODUCER = as.character(PRODUCER),
         Same_class = (PRODUCER == Producer_no2020))

# Summary: how many trees changed class?
cat("Producer classification with vs. without 2020:\n")
```

    ## Producer classification with vs. without 2020:

``` r
cat("  Total trees:", nrow(prod_compare), "\n")
```

    ##   Total trees: 39

``` r
cat("  Same class:", sum(prod_compare$Same_class, na.rm = TRUE), "\n")
```

    ##   Same class: 9

``` r
cat("  Changed class:", sum(!prod_compare$Same_class, na.rm = TRUE), "\n\n")
```

    ##   Changed class: 30

``` r
# Cross-tabulation
cross_tab <- table(Original = prod_compare$PRODUCER,
                   No2020   = prod_compare$Producer_no2020)
cat("Cross-tabulation (rows = original, cols = without 2020):\n")
```

    ## Cross-tabulation (rows = original, cols = without 2020):

``` r
print(cross_tab)
```

    ##         No2020
    ## Original Super Medium Low
    ##     High     6      0   0
    ##     Low      0      3   9
    ##     Mid      8     10   3

``` r
# Trees that changed
changed <- prod_compare %>%
  filter(!Same_class) %>%
  select(Tree, PRODUCER, Producer_no2020, mean_count_no2020)
if (nrow(changed) > 0) {
  cat("\nTrees that changed classification:\n")
  print(changed)
} else {
  cat("\n✓ All trees retain the same producer classification when 2020 is excluded.\n")
}
```

    ## 
    ## Trees that changed classification:
    ## # A tibble: 30 × 4
    ##    Tree  PRODUCER Producer_no2020 mean_count_no2020
    ##    <chr> <chr>    <fct>                       <dbl>
    ##  1 1     High     Super                      141   
    ##  2 2     Mid      Medium                      27.7 
    ##  3 3     Mid      Low                          7.33
    ##  4 4     Mid      Medium                      10.7 
    ##  5 8     Mid      Super                       43   
    ##  6 11    Mid      Super                       36.3 
    ##  7 12    Mid      Super                       34.3 
    ##  8 13    High     Super                      132   
    ##  9 14    Mid      Medium                      19.7 
    ## 10 15    High     Super                      127   
    ## # ℹ 20 more rows

``` r
knitr::kable(prod_compare %>%
               select(Tree, Original = PRODUCER, No2020 = Producer_no2020,
                      Mean_no2020 = mean_count_no2020, Same_class) %>%
               arrange(Original, No2020),
             caption = "**Producer classification with and without 2020.**
             Original: classification from wytham_full_dataset.xlsx (all years).
             No2020: reclassified on 2021–2023 means using tertile thresholds.
             Same_class = TRUE if the producer category did not change.",
             digits = 1)
```

| Tree | Original | No2020 | Mean_no2020 | Same_class |
|:-----|:---------|:-------|------------:|:-----------|
| 1    | High     | Super  |       141.0 | FALSE      |
| 13   | High     | Super  |       132.0 | FALSE      |
| 15   | High     | Super  |       127.0 | FALSE      |
| 33   | High     | Super  |       178.3 | FALSE      |
| 36   | High     | Super  |       212.0 | FALSE      |
| 39   | High     | Super  |       136.0 | FALSE      |
| 25   | Low      | Medium |        14.0 | FALSE      |
| 32   | Low      | Medium |         8.0 | FALSE      |
| 37   | Low      | Medium |        11.0 | FALSE      |
| 5    | Low      | Low    |         2.7 | TRUE       |
| 6    | Low      | Low    |         5.3 | TRUE       |
| 7    | Low      | Low    |         0.0 | TRUE       |
| 9    | Low      | Low    |         1.3 | TRUE       |
| 16   | Low      | Low    |         2.0 | TRUE       |
| 17   | Low      | Low    |         1.3 | TRUE       |
| 21   | Low      | Low    |         0.0 | TRUE       |
| 27   | Low      | Low    |         1.7 | TRUE       |
| 38   | Low      | Low    |         5.0 | TRUE       |
| 8    | Mid      | Super  |        43.0 | FALSE      |
| 11   | Mid      | Super  |        36.3 | FALSE      |
| 12   | Mid      | Super  |        34.3 | FALSE      |
| 19   | Mid      | Super  |        57.0 | FALSE      |
| 20   | Mid      | Super  |        60.0 | FALSE      |
| 23   | Mid      | Super  |        40.7 | FALSE      |
| 30   | Mid      | Super  |        49.3 | FALSE      |
| 35   | Mid      | Super  |        34.3 | FALSE      |
| 2    | Mid      | Medium |        27.7 | FALSE      |
| 4    | Mid      | Medium |        10.7 | FALSE      |
| 14   | Mid      | Medium |        19.7 | FALSE      |
| 22   | Mid      | Medium |        15.7 | FALSE      |
| 24   | Mid      | Medium |         9.7 | FALSE      |
| 26   | Mid      | Medium |        13.0 | FALSE      |
| 28   | Mid      | Medium |        17.7 | FALSE      |
| 31   | Mid      | Medium |        28.0 | FALSE      |
| 34   | Mid      | Medium |        32.3 | FALSE      |
| 41   | Mid      | Medium |        29.7 | FALSE      |
| 3    | Mid      | Low    |         7.3 | FALSE      |
| 18   | Mid      | Low    |         3.7 | FALSE      |
| 29   | Mid      | Low    |         7.3 | FALSE      |

**Producer classification with and without 2020.** Original:
classification from wytham_full_dataset.xlsx (all years). No2020:
reclassified on 2021–2023 means using tertile thresholds. Same_class =
TRUE if the producer category did not change.

``` r
write.csv(prod_compare %>%
            select(Tree, TREE_ID, Original_PRODUCER = PRODUCER,
                   Producer_no2020, mean_count_no2020, sum_count_no2020,
                   Same_class),
          file.path(TAB_OUT, "Superproducer_without_2020.csv"),
          row.names = FALSE)
```

------------------------------------------------------------------------

# Session information

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26200)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.utf8 
    ## [2] LC_CTYPE=English_United Kingdom.utf8   
    ## [3] LC_MONETARY=English_United Kingdom.utf8
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.utf8    
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.3.2    scales_1.4.0       car_3.1-5          carData_3.0-6     
    ##  [5] plotrix_3.8-14     geomtextpath_0.2.0 viridis_0.6.5      viridisLite_0.4.2 
    ##  [9] ggthemes_5.2.0     readxl_1.4.5       lubridate_1.9.4    forcats_1.0.1     
    ## [13] stringr_1.6.0      dplyr_1.1.4        purrr_1.2.2        readr_2.1.6       
    ## [17] tidyr_1.3.1        tibble_3.3.0       ggplot2_4.0.1      tidyverse_2.0.0   
    ## [21] rmarkdown_2.30    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.6         generics_0.1.4     lattice_0.22-7     stringi_1.8.7     
    ##  [5] hms_1.1.4          digest_0.6.39      magrittr_2.0.4     evaluate_1.0.5    
    ##  [9] grid_4.5.2         timechange_0.3.0   RColorBrewer_1.1-3 fastmap_1.2.0     
    ## [13] Matrix_1.7-4       cellranger_1.1.0   Formula_1.2-5      gridExtra_2.3     
    ## [17] mgcv_1.9-3         textshaping_1.0.5  abind_1.4-8        cli_3.6.5         
    ## [21] rlang_1.2.0        splines_4.5.2      withr_3.0.2        yaml_2.3.12       
    ## [25] otel_0.2.0         tools_4.5.2        tzdb_0.5.0         vctrs_0.7.3       
    ## [29] R6_2.6.1           lifecycle_1.0.5    ragg_1.5.1         pkgconfig_2.0.3   
    ## [33] pillar_1.11.1      gtable_0.3.6       glue_1.8.0         systemfonts_1.3.2 
    ## [37] xfun_0.57          tidyselect_1.2.1   knitr_1.51         farver_2.1.2      
    ## [41] nlme_3.1-168       htmltools_0.5.9    labeling_0.4.3     compiler_4.5.2    
    ## [45] S7_0.2.1
