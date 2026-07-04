# QPM-GWAS

**Genome-Wide Association Studies for Kernel Size and Number in Quality Protein Maize (*Zea mays* L.)**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-4.3.3-blue.svg)](https://cran.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1007%2Fs13353--026--01078--3-green.svg)](https://doi.org/10.1007/s13353-026-01078-3)

---

## Background

Kernel size and number are primary determinants of yield in maize, yet their genetic architecture in **Quality Protein Maize (QPM)** — a biofortified variety critical for food and nutritional security in South Asia and Sub-Saharan Africa — remains poorly characterised. QPM delivers nearly twice the usable protein of conventional maize due to the *opaque2* mutation, but yield gaps continue to limit its adoption by smallholder farmers.

This repository documents the complete analytical pipeline used to identify the genomic regions governing eight yield-related traits in a diverse panel of **149 QPM inbred lines** evaluated across multiple agroclimatic environments in India. The study successfully identified **27 stable marker-trait associations (MTAs)**, including **9 pleiotropic loci** and novel genomic bins on chromosomes 3 and 5, providing a foundation for genomics-assisted yield improvement in QPM.

---

## Published Paper

> **Rainy R.**, Sandhu S.S., Kumar S., Kumar R., Vikal Y., Sharma P. (2026).  
> Genome-wide Association in Quality Protein Maize: Identifying Key Loci for Yield Enhancement.  
> *Journal of Applied Genetics.*  
> DOI: [10.1007/s13353-026-01078-3](https://doi.org/10.1007/s13353-026-01078-3)

---

## Pipeline Overview

The pipeline follows a structured, sequential workflow from raw field data to candidate gene annotation:

```
Raw Phenotypic Data
        │
        ▼
00_phenotypic_analysis.R
   Descriptive stats · ANOVA (lattice) · Heritability · BLUEs · Correlation
        │
        │── GxE ANOVA (SPSS Statistics)
        │   Multi-environment interaction analysis
        │
        ▼
01_QC.R
   SNP filtering: MAF (5%) · Heterozygosity (20%) · Missingness (20%)
   Tools: TASSEL · PLINK
        │
        ▼
02_imputation.sh
   Missing genotype imputation
   Tool: Beagle v5.5
        │
        ▼
03_GWAS.R
   Association mapping: FarmCPU · MLM
   Package: GAPIT3 (R)
        │
        ▼
04_visualization.R
   Manhattan plots · QQ plots · LD decay curves
        │
        ▼
05_candidate_genes.R
   LD interval extraction · Gene annotation
   Databases: Ensembl Plants · MaizeGDB · STRING · InterPro · DAVID
```

---

## Repository Structure

```
QPM-GWAS/
│
├── README.md
│
├── scripts/
│   ├── 00_phenotypic_analysis.R   # ✅ Available — full phenotypic pipeline
│   ├── 01_QC.R                    # 🔄 In progress
│   ├── 02_imputation.sh           # 🔄 In progress
│   ├── 03_GWAS.R                  # 🔄 In progress
│   ├── 04_visualization.R         # 🔄 In progress
│   └── 05_candidate_genes.R       # 🔄 In progress
│
├── data/
│   └── README.md                  # Data description and access notes
│
├── figures/
│   └── README.md                  # Description of output figures
│
└── results/
    └── README.md                  # Description of output files
```

---

## Script 00 — Phenotypic Analysis

### What it does

A generalised, reusable pipeline for phenotypic analysis of multi-trait data from plant breeding trials conducted in a **resolvable incomplete block (lattice) design**. Although developed for QPM, the script works with any crop and any set of traits with minimal configuration.

| Section | Analysis |
|---------|----------|
| 1 | Package setup and data import (Excel or CSV) |
| 2 | Data validation and trait auto-detection |
| 3 | Descriptive statistics — mean, SD, CV%, skewness, kurtosis |
| 4 | Compact scattered boxplots per trait |
| 5 | Frequency distributions with mean and median lines |
| 6 | ANOVA for lattice design with outlier detection |
| 7 | Critical Difference (CD) at 5% |
| 8 | Genetic variability — GV, PV, h², GA, GA% |
| 9 | Normality testing, Box-Cox transformation, BLUE estimation |
| 10 | Pearson correlation — corrplot and ggpairs panel |

### Experimental design

```
149 QPM inbred lines  ×  2 replications  ×  multiple environments
Resolvable incomplete block (lattice) design
Model: trait ~ Genotype + Replication + Replication:Block

GxE interaction analysis performed separately in SPSS Statistics
```

### Traits analysed

| Code | Trait | Unit |
|------|-------|------|
| KL | Kernel Length | mm |
| KW | Kernel Width | mm |
| KT | Kernel Thickness | mm |
| CD | Cob Diameter | mm |
| CL | Cob Length | cm |
| KNPR | Kernel Number per Row | count |
| RN | Row Number | count |
| KNPC | Kernel Number per Cob | count |

### Key outputs

| File | Contents |
|------|----------|
| `results/Descriptive_Stats.xlsx` | Summary statistics and CV table |
| `results/ANOVA_CD.xlsx` | ANOVA table and Critical Difference values |
| `results/Genetic_Variability.xlsx` | GV, PV, h², GA, GA% per trait |
| `results/BLUEs_Normality.xlsx` | BLUEs for GWAS input + normality check |
| `results/Correlation.xlsx` | Pearson correlation matrix and p-values |
| `figures/Boxplot_<trait>.png` | Compact scattered boxplot per trait |
| `figures/FreqDist_<trait>.png` | Frequency distribution per trait |
| `figures/Genetic_Variability.png` | h² and GA% bar plots |
| `figures/Correlation_corrplot.png` | Upper-triangle correlation heatmap |
| `figures/Correlation_ggpairs.tiff` | Full pairwise correlation panel |

---

## How to Use

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/QPM-GWAS.git
cd QPM-GWAS
```

### 2. Prepare your data

Format your input file as one row per plot with columns for:
- Genotype identifier
- Replication
- Incomplete block
- Numeric trait values

Both `.xlsx` and `.csv` formats are supported.

### 3. Configure Section 0

Open `scripts/00_phenotypic_analysis.R` and edit only the **USER CONFIGURATION** block at the top:

```r
INPUT_FILE      <- "data/your_data.xlsx"   # path to your file
SHEET           <- 1                        # Excel sheet number (ignored for CSV)
COL_GENOTYPE    <- "Genotype"              # your genotype column name
COL_REPLICATION <- "Rep"                   # your replication column name
COL_BLOCK       <- "Block"                 # your block column name
TRAITS_MANUAL   <- NULL                    # NULL = auto-detect traits
N_REPS          <- 2                       # number of replications
```

### 4. Run the script

```r
source("scripts/00_phenotypic_analysis.R")
```

All outputs are written automatically to `results/` and `figures/`.

---

## Requirements

- **R** >= 4.3.3
- The script installs all required packages automatically on first run.

Key packages: `lme4`, `emmeans`, `MASS`, `ggplot2`, `GGally`, `corrplot`,
`cowplot`, `Hmisc`, `openxlsx`, `readxl`, `moments`, `psych`

---

## Data Availability

Raw phenotypic and genotypic data are available upon reasonable request to
the corresponding author:

**Dr. Priti Sharma**  
Assistant Professor, School of Agricultural Biotechnology  
Punjab Agricultural University, Ludhiana, India  
📧 pritisharma@pau.edu

---

## Key Findings

| Parameter | Value |
|-----------|-------|
| Panel size | 149 QPM inbred lines |
| Environments | Multiple agroclimatic locations, India |
| Stable MTAs identified | 27 |
| Pleiotropic loci | 9 |
| PVE range | 5.64 – 27.99% |
| Novel genomic bins | Chr 3 (3.03–3.04) and Chr 5 (5.01) |

---

## Citation

If you use this pipeline or adapt the scripts for your own research, please cite:

```
Rainy R., Sandhu S.S., Kumar S., Kumar R., Vikal Y., Sharma P. (2026).
Genome-wide Association in Quality Protein Maize: Identifying Key Loci
for Yield Enhancement. Journal of Applied Genetics.
DOI: 10.1007/s13353-026-01078-3
```

---

## Contact

**Rainy** — rainy122001@gmail.com · [LinkedIn](https://www.linkedin.com/in/rainy)

---

## License

This project is licensed under the MIT License.  
See the [LICENSE](LICENSE) file for details.
