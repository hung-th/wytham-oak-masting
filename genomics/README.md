# Genomic data processing (computing cluster)

The steps below were run on a Linux computing cluster. Their outputs are the
genomic inputs stored in `data/raw/genomics/`, which the notebook
`analysis/04_genomic_associations.Rmd` uses.

## 1. Reads to filtered SNPs

Whole-genome Illumina NovaSeq 6000 reads (150 bp paired-end) for 37 trees
were processed as follows:

| Step | Tool | Settings |
|---|---|---|
| Adapter and quality trimming | fastp 0.24.0 | |
| Mapping to *Q. robur* PM1N (Plomion et al. 2018) | bwa-mem2 2.2.1 | |
| Sorting | samtools | |
| Pileup and variant calling | bcftools 1.21 | |
| SNP filtering | VCFtools 0.1.16 | MAF ≥ 0.01, minimum depth 50, < 10% missing data |
| Indel-proximal SNPs | | SNPs within 150 bp of an indel removed |
| LD pruning | PLINK 1.90b6.21 | `--indep-pairwise 50 5 0.2` (50-SNP window, 5-SNP step, r² 0.2) |

Read and mapping statistics per tree are in
`data/raw/genomics/sequencing_statistics.csv` (Supplementary Table 12).

## 2. Population structure and kinship

| Output (`data/raw/genomics/`) | Tool | Input |
|---|---|---|
| `population_structure/pca.eigenval`, `pca.eigenvec` | PLINK `--pca` | LD-pruned SNPs |
| `population_structure/kinship_full_snpset.relatedness2` | VCFtools `--relatedness2` (KING-robust φ) | filtered SNPs before LD pruning |
| `genotypes/WWO_37WGS_pruned.012`, `.012.indv`, `.012.pos` | VCFtools `--012` | LD-pruned SNPs |
| `genotypes/WWO_37WGS_pruned.lfmm` | 012 matrix in LFMM format (missing calls = 9) | LD-pruned SNPs |

## 3. Scripts in this folder

| Script | Purpose | Output |
|---|---|---|
| `run_snmf.R` | sNMF (LEA) on 100,000 randomly sampled LD-pruned SNPs, K = 1–6, three repetitions | `snmf_Q.rds` (K = 1), cross-entropies (`cache/gwas/snmf_crossentropy.rds`) |
| `gwas_impute_parallel.R` | Imputes missing genotypes by regressing each SNP on the sNMF ancestry matrix | `genotypes/WWO_37WGS_pruned_imputed.lfmm`, `.012.pos` |
| `gwas_lfmm.R` | LFMM2 association tests (K = 1) with genome-wide genomic control and Q-values | `cache/gwas/lfmm2_results.rds` |

`gwas_lfmm.R` runs in a few minutes on a desktop computer and is called
automatically by `analysis/04_genomic_associations.Rmd` when its cached result
is missing or out of date.
