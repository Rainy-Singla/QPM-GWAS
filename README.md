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
   SNP filtering: MAF (5%) · Missingness (20%) · Heterozygosity (20%, TASSEL)
   Tools: PLINK · TASSEL
        │
        ▼
02_imputation.sh
   Missing genotype imputation
   Tool: Beagle v5.5
        │
        ▼
03_GWAS.R
   LD decay · Marker density · Association mapping
   Models: FarmCPU · MLM · BLINK · CMLM · Super
   Package: GAPIT3 (R)
        │
        ▼
04_visualization.R
   Manhattan plot · Venn diagram
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
│   ├── 01_QC.R                    # ✅ Available — SNP quality control
│   ├── 02_imputation.sh           # ✅ Available — Beagle v5.5 imputation
│   ├── 03_GWAS.R                  # ✅ Available — LD decay, marker density, GAPIT3 GWAS
│   ├── 04_visualization.R         # ✅ Available — Manhattan plot, Venn diagram
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
 
## Scripts

### 00 — Phenotypic Analysis

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

**Experimental design:**
```
149 QPM inbred lines  ×  2 replications  ×  multiple environments
Resolvable incomplete block (lattice) design
Model: trait ~ Genotype + Replication + Replication:Block
GxE interaction analysis performed separately in SPSS Statistics
```

**Traits analysed:**

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

---

### 01 — SNP Quality Control

Performs QC filtering of GBS-derived SNP data using PLINK. Heterozygosity filtering (20%) was performed upstream in TASSEL prior to running this script.

| Step | Filter | Threshold |
|------|--------|-----------|
| 1 | Heterozygosity (TASSEL) | ≤ 20% per SNP |
| 2 | Minor Allele Frequency | ≥ 5% |
| 3 | SNP missingness | ≤ 20% |

**Prerequisites:** PLINK v1.9 — [download here](https://www.cog-genomics.org/plink/)

---

### 02 — Imputation

Missing genotype imputation using Beagle v5.5.

| Parameter | Value |
|-----------|-------|
| Tool | Beagle v5.5 |
| Memory allocation | 8GB (-Xmx8g) |
| Input | QC-filtered VCF from 01_QC.R |
| Output | Imputed VCF (.vcf.gz), decompressed .vcf |

**Prerequisites:** Java >= 8 and Beagle v5.5 JAR — [download here](https://faculty.washington.edu/browning/beagle/beagle.html)

---

### 03 — GWAS
 
LD decay estimation, chromosome marker density plot, and association mapping via GAPIT3.
 
| Section | Analysis |
|---------|----------|
| 1 | LD decay curve — Hill-Weir model fitted via NLS (Remington et al., 2001) |
| 2 | Chromosome marker density plot (CMplot, 1 Mb bins) |
| 3 | GWAS — FarmCPU, MLM, BLINK, CMLM, Super via GAPIT3 |
| 4 | Significant SNP filtering across models at -log10(p) threshold |
 
**NOTE — Model selection:** All five models run by default. Edit `GWAS_MODELS` in Section 0 to run specific models only. For stable MTA identification, SNPs consistently detected by FarmCPU and MLM are recommended as a conservative approach.
 
**Prerequisites:** GAPIT3 installed from GitHub — `devtools::install_github("jiabowang/GAPIT3")`
 
---
 
### 04 — Visualisation
 
Manhattan plot and Venn diagram of significant SNPs across GWAS models.
 
| Section | Analysis |
|---------|----------|
| 1 | Manhattan plot — genome-wide p-value distribution with significant SNPs highlighted (CMplot) |
| 2 | Venn diagram — overlap of significant SNPs across models (ggvenn) |
 
**NOTE:** An UpSet plot alternative (UpSetR) is included as commented code for datasets with more than four models or complex overlaps.
 
---
 
### 05 — Candidate Gene Annotation *(in progress)*

LD interval extraction and functional annotation via Ensembl Plants, MaizeGDB, STRING, InterPro, and DAVID.

---

## How to Use

### 1. Clone the repository

```bash
git clone https://github.com/Rainy-Singla/QPM-GWAS.git
cd QPM-GWAS
```

### 2. Prepare your data

Format your input file as one row per plot with columns for:
- Genotype identifier
- Replication
- Incomplete block
- Numeric trait values

Both `.xlsx` and `.csv` formats are supported for phenotypic analysis.

### 3. Configure Section 0

Each script has a **USER CONFIGURATION** block at the top — edit only this section to point to your files and set your parameters. No other changes are needed.

### 4. Run scripts in order

```r
source("scripts/00_phenotypic_analysis.R")
source("scripts/01_QC.R")
# bash scripts/02_imputation.sh
source("scripts/03_GWAS.R")
source("scripts/04_visualization.R")
```

All outputs are written automatically to `results/` and `figures/`.

---

## Requirements

- **R** >= 4.3.3
- **PLINK** v1.9 (for `01_QC.R`)
- **Java** >= 8 and **Beagle** v5.5 JAR (for `02_imputation.sh`)
- **TASSEL** >= 5.0 (for heterozygosity filtering)

All R packages are installed automatically on first run.

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

## Data Availability

Raw phenotypic and genotypic data are available upon reasonable request to
the corresponding author:

**Dr. Priti Sharma**  
Assistant Professor, School of Agricultural Biotechnology  
Punjab Agricultural University, Ludhiana, India  
📧 pritisharma@pau.edu

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
