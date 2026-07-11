# =============================================================================
# GWAS Pipeline: LD Decay, Marker Density, and Association Mapping
# =============================================================================
#
# Description:
#   This script performs genome-wide association analysis on imputed SNP data.
#   The pipeline covers:
#
#     1. LD decay estimation and visualisation
#        (Hill-Weir model; Remington et al., 2001)
#     2. Chromosome marker density plot (CMplot)
#     3. GWAS using GAPIT3 — five statistical models:
#          FarmCPU, MLM, BLINK, CMLM, Super
#     4. Significant SNP filtering across models
#
# NOTE — Model selection:
#   All five models are run by default. If you want to run only specific
#   models, edit the GWAS_MODELS vector in Section 0 and keep only the
#   models you need. For example:
#     GWAS_MODELS <- c("FarmCPU", "MLM")   # run only FarmCPU and MLM
#   Recommended for initial screening: FarmCPU + MLM
#   Recommended for stable MTA identification: SNPs significant in both
#   FarmCPU and MLM (conservative approach used in the published study)
#
# Input:
#   - TASSEL LD output file (.txt)       — for LD decay
#   - QC-filtered VCF (.vcf)             — for marker density plot
#   - Imputed HapMap file (.hmp.txt)     — genotype input for GAPIT
#   - Phenotype file (.txt)              — BLUE/BLUP values per genotype
#
# Output:
#   - figures/LD_decay.pdf               — LD decay curve
#   - figures/Marker_Density.jpg         — chromosome marker density
#   - results/GAPIT_<model>_<trait>/     — GAPIT output folders per model
#   - results/significant_SNPs_<model>.xlsx — filtered significant SNPs
#
# References:
#   Remington et al. (2001) PNAS 98(20):11479-11484
#   Wang et al. (2021) Genomics, Proteomics & Bioinformatics (GAPIT3)
#
# Author:  Rainy
# Contact: rainy122001@gmail.com
# DOI:     https://doi.org/10.1007/s13353-026-01078-3
# License: MIT
# =============================================================================


# =============================================================================
# 0. USER CONFIGURATION — edit this section before running
# =============================================================================

# --- 0.1 Input files ---------------------------------------------------------

# TASSEL LD output file (tab-delimited, with R^2 and distance columns)
LD_FILE <- "data/LD_output.txt"

# QC-filtered VCF file (from 01_QC.R, before imputation)
VCF_FILE <- "data/final_filtered_vcf.vcf"

# Imputed HapMap genotype file (output from Beagle, converted to HapMap format)
HAPMAP_FILE <- "data/filtered_imputed.hmp.txt"

# Phenotype file: tab-delimited, first column = genotype ID (taxa),
# remaining columns = trait BLUEs/BLUPs (one column per trait)
PHENOTYPE_FILE <- "data/BLUEs.txt"

# Trait name — used for labelling output files
# Should match the column name in your phenotype file
TRAIT_NAME <- "KNPC"

# --- 0.2 GWAS models ---------------------------------------------------------
# All five models run by default. Remove any you do not need.
# Available: "FarmCPU", "MLM", "BLINK", "CMLM", "Super"

GWAS_MODELS <- c("FarmCPU", "MLM", "BLINK", "CMLM", "Super")

# Number of PCs to include for population structure correction
PCA_TOTAL <- 3

# --- 0.3 LD decay parameters -------------------------------------------------

# Number of genotypes used in the LD analysis (N in Hill-Weir formula)
# This should match the number of individuals in your dataset
N_GENOTYPES <- 149

# LD decay half-decay threshold — the r² value at which to mark decay
LD_HALF_DECAY_THRESHOLD <- 0.1

# --- 0.4 Significant SNP filtering -------------------------------------------

# -log10(p) threshold for significance
# 3.5 corresponds to p < ~3.2e-4 (commonly used in plant GWAS)
# Adjust based on your sample size and multiple testing approach
LOG10P_THRESHOLD <- 3.5

# --- 0.5 Output paths --------------------------------------------------------

OUT_LD_PLOT      <- "figures/LD_decay.pdf"
OUT_GWAS_RESULTS <- "results/"

# Set working directory — update to your local project folder
# setwd("path/to/your/project")


# =============================================================================
# 1. SETUP
# =============================================================================

required_packages <- c(
  "ggplot2", "dplyr", "data.table", "openxlsx",
  "readxl", "CMplot"
)

new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages) > 0) {
  message("Installing missing packages: ", paste(new_packages, collapse = ", "))
  install.packages(new_packages)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(openxlsx)
  library(readxl)
  library(CMplot)
})

# Install GAPIT3 from GitHub if not already installed
if (!requireNamespace("GAPIT", quietly = TRUE)) {
  message("Installing GAPIT3 from GitHub...")
  if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
  devtools::install_github("jiabowang/GAPIT3", force = TRUE)
}
library(GAPIT)

# Create output directories
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("=== GWAS Pipeline ===\n")
cat("Trait:", TRAIT_NAME, "\n")
cat("Models:", paste(GWAS_MODELS, collapse = ", "), "\n\n")


# =============================================================================
# 2. LD DECAY
# =============================================================================
# Fits the Hill-Weir (1988) expectation curve to observed r² vs. distance
# data from TASSEL. The non-linear least squares (NLS) model estimates
# the recombination parameter rho (4Nr), which is used to generate a
# smooth expected LD decay curve.
#
# Reference:
#   Remington et al. (2001) PNAS 98(20):11479-11484
#   https://doi.org/10.1073/pnas.201394398
#
# Input: TASSEL LD output — calculate with sliding window (e.g. 50 SNPs)
# at MAF >= 0.05 in TASSEL's LD analysis module.

cat("--- Section 2: LD Decay ---\n")

# Import TASSEL LD output
ld <- read.delim(LD_FILE, stringsAsFactors = FALSE, header = TRUE, sep = "\t")

# Remove rows with NaN in R² or distance
ld_clean <- ld[ld$R.2 != "NaN", ]
ld_clean$dist <- as.numeric(ld_clean$Dist_bp)
ld_clean <- ld_clean[!is.na(ld_clean$dist), ]
ld_clean$rsq <- as.numeric(ld_clean$R.2)

# Select relevant columns: Locus1, Locus2, dist, rsq, and N
# Adjust column indices if your TASSEL output has a different structure
ld_file <- ld_clean[, c(1, 2, 7, 8, 15:19)]
ld_file$N <- N_GENOTYPES   # set N to number of genotypes in your panel

# Fit Hill-Weir NLS model
# C (rho = 4Nr) is the recombination parameter; start value of 0.1 works
# for most datasets — adjust if model fails to converge
Cstart <- c(C = 0.1)

modelC <- nls(
  rsq ~ ((10 + C * dist) / ((2 + C * dist) * (11 + C * dist))) *
    (1 + ((3 + C * dist) * (12 + 12 * C * dist + (C * dist)^2)) /
       (2 * N * (2 + C * dist) * (11 + C * dist))),
  data    = ld_file,
  start   = Cstart,
  control = nls.control(maxiter = 100)
)

# Extract fitted rho
rho <- summary(modelC)$parameters[1]
cat("  Fitted rho (4Nr):", round(rho, 6), "\n")

# Generate expected LD decay curve
newrsq <- ((10 + rho * ld_file$dist) / ((2 + rho * ld_file$dist) *
             (11 + rho * ld_file$dist))) *
  (1 + ((3 + rho * ld_file$dist) *
          (12 + 12 * rho * ld_file$dist + (rho * ld_file$dist)^2)) /
     (2 * ld_file$N * (2 + rho * ld_file$dist) *
        (11 + rho * ld_file$dist)))

newfile <- data.frame(dist = ld_file$dist, newrsq = newrsq)
newfile <- newfile[order(newfile$dist), ]

# Calculate half-decay distance at LD_HALF_DECAY_THRESHOLD
maxld        <- max(newfile$newrsq, na.rm = TRUE)
halfdecay    <- LD_HALF_DECAY_THRESHOLD
halfdecaydist <- newfile$dist[which.min(abs(newfile$newrsq - halfdecay))]
cat("  Half-decay distance at r² =", LD_HALF_DECAY_THRESHOLD,
    ":", round(halfdecaydist, 0), "bp\n")

# Plot LD decay
pdf(OUT_LD_PLOT, height = 5, width = 5)
par(mar = c(5, 8, 4, 2) + 0.1)
plot(ld_file$dist, ld_file$rsq,
     pch  = ".", cex = 2,
     xlab = "Distance (bp)",
     ylab = expression(LD ~ (r^2)),
     col  = "grey")
lines(newfile$dist, newfile$newrsq, col = "red", lwd = 2)
abline(h = LD_HALF_DECAY_THRESHOLD, col = "blue")
abline(v = halfdecaydist, col = "green")
mtext(round(halfdecaydist, 0), side = 1, line = 0.05,
      at = halfdecaydist, cex = 0.75, col = "green")
dev.off()
cat("  Saved:", OUT_LD_PLOT, "\n\n")


# =============================================================================
# 3. CHROMOSOME MARKER DENSITY PLOT
# =============================================================================
# Visualises the distribution of SNPs across chromosomes in 1 Mb windows.
# Uses the QC-filtered VCF (before imputation) to show post-QC SNP density.

cat("--- Section 3: Chromosome Marker Density ---\n")

# Read VCF, skipping header lines starting with "##"
vcf_data <- read.table(VCF_FILE, comment.char = "#",
                        header = FALSE, sep = "\t")

# Standard VCF columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO...
colnames(vcf_data)[1:3] <- c("Chromosome", "Position", "SNP")

# Keep only SNP, chromosome, and position
snp_info <- vcf_data[, c("SNP", "Chromosome", "Position")]

# Remove "chr" prefix from chromosome names for CMplot compatibility
snp_info$Chromosome <- gsub("chr", "", snp_info$Chromosome, ignore.case = TRUE)
snp_info$Chromosome <- as.numeric(as.character(snp_info$Chromosome))

# Remove SNPs on unassigned scaffolds (NA chromosome after conversion)
snp_info <- snp_info[!is.na(snp_info$Chromosome), ]

cat("  SNPs for density plot:", nrow(snp_info), "\n")

CMplot(snp_info,
       plot.type   = "d",
       bin.size    = 1e6,
       chr.den.col = c("gray", "green", "yellow", "orange", "red"),
       file        = "jpg",
       dpi         = 300,
       file.output = TRUE,
       output.file  = "figures/Marker_Density",
       width       = 10,
       height      = 6)

cat("  Saved: figures/Marker_Density.jpg\n\n")


# =============================================================================
# 4. GWAS — GAPIT3
# =============================================================================
# Runs genome-wide association mapping using GAPIT3.
#
# Models:
#   FarmCPU — Fixed and random model Circulating Probability Unification
#              Controls false positives while maintaining power
#   MLM     — Mixed Linear Model; uses kinship + PCA for correction
#   BLINK   — Bayesian information and Linkage disequilibrium Iteratively
#              Nested Keyway; fast and powerful for large datasets
#   CMLM    — Compressed MLM; groups similar individuals for efficiency
#   Super   — SuperBLUP; optimal for highly structured populations
#
# For stable MTA identification, report SNPs consistently detected by
# two or more models (e.g. FarmCPU + MLM as in the published study).
# Adjust based on your population structure and research objective.
#
# Input formats:
#   Phenotype: tab-delimited .txt, first column = taxa (genotype IDs),
#              remaining columns = trait values (one per trait)
#   Genotype:  HapMap format (.hmp.txt) — convert from VCF using TASSEL:
#              File > Load Data > VCF, then Data > Export > Hapmap

cat("--- Section 4: GWAS (GAPIT3) ---\n")
cat("  Models:", paste(GWAS_MODELS, collapse = ", "), "\n")

# Import phenotype and genotype data
myY <- read.delim(PHENOTYPE_FILE, head = TRUE)
myG <- read.delim(HAPMAP_FILE,    head = FALSE)

cat("  Genotypes loaded:", nrow(myY), "\n")
cat("  SNPs loaded:", nrow(myG) - 1, "\n")
cat("  Running GAPIT... (this may take several minutes)\n")

# Create output folder for GAPIT results
gapit_out <- file.path(OUT_GWAS_RESULTS, paste0("GAPIT_", TRAIT_NAME))
dir.create(gapit_out, showWarnings = FALSE, recursive = TRUE)
setwd(gapit_out)

# Run GAPIT with all specified models
myGAPIT <- GAPIT(
  Y         = myY,
  G         = myG,
  PCA.total = PCA_TOTAL,
  model     = GWAS_MODELS
)

# Return to project root after GAPIT (GAPIT changes working directory)
setwd("../../")
cat("  GAPIT complete. Results saved to:", gapit_out, "\n\n")


# =============================================================================
# 5. SIGNIFICANT SNP FILTERING
# =============================================================================
# Filters GAPIT output for each model at the specified -log10(p) threshold.
# Saves significant SNPs per model to individual Excel files.
#
# NOTE: GAPIT names output files as:
#   GAPIT.Association.GWAS_Results.<Model>.<Trait>.csv
# Adjust the file_paths pattern below if your GAPIT version differs.

cat("--- Section 5: Significant SNP Filtering ---\n")
cat("  Threshold: -log10(p) >=", LOG10P_THRESHOLD, "\n")

sig_snps_all <- data.frame()   # collects significant SNPs across all models

for (model in GWAS_MODELS) {
  # Construct expected GAPIT output filename
  file_name <- paste0(
    gapit_out, "/GAPIT.Association.GWAS_Results.",
    model, ".", TRAIT_NAME, ".csv"
  )

  if (!file.exists(file_name)) {
    warning("File not found, skipping: ", file_name)
    next
  }

  gwas_data <- fread(file_name)

  # Filter significant SNPs
  sig <- gwas_data %>%
    filter(!is.na(P.value)) %>%
    mutate(
      logP  = -log10(as.numeric(P.value)),
      Model = model
    ) %>%
    filter(logP >= LOG10P_THRESHOLD)

  cat("  ", model, "— significant SNPs:", nrow(sig), "\n")

  # Save per-model results
  out_file <- file.path(
    OUT_GWAS_RESULTS,
    paste0("significant_SNPs_", model, "_", TRAIT_NAME, ".xlsx")
  )
  write.xlsx(sig, out_file)
  cat("  Saved:", out_file, "\n")

  sig_snps_all <- rbind(sig_snps_all, sig)
}

# Save combined significant SNPs across all models
if (nrow(sig_snps_all) > 0) {
  combined_out <- file.path(
    OUT_GWAS_RESULTS,
    paste0("significant_SNPs_ALL_models_", TRAIT_NAME, ".xlsx")
  )
  write.xlsx(sig_snps_all, combined_out)
  cat("\nCombined significant SNPs saved to:", combined_out, "\n")
  cat("Total significant SNPs (with duplicates across models):",
      nrow(sig_snps_all), "\n")
  cat("Unique significant SNPs:",
      length(unique(sig_snps_all$SNP)), "\n")
}

cat("\n=== GWAS pipeline complete. Proceed to 04_visualization.R ===\n")
sessionInfo()
