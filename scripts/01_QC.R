# =============================================================================
# SNP Quality Control (QC) Pipeline
# =============================================================================
#
# Description:
#   This script performs quality control on raw GBS-derived genotypic data
#   in VCF format using PLINK. The pipeline converts the VCF to PLINK binary
#   format, applies QC filters, and exports a filtered VCF for downstream
#   imputation and association mapping.
#
#   QC filters applied:
#     - Minor Allele Frequency (MAF)    : removes rare variants
#     - SNP missingness (--geno)        : removes poorly genotyped SNPs
#     - Heterozygosity                  : performed externally in TASSEL
#                                         (see note below)
#
# Prerequisites:
#   PLINK v1.9 must be installed and accessible from the command line.
#   Download: https://www.cog-genomics.org/plink/
#   Place plink.exe (Windows) or plink (Linux/Mac) in your working directory
#   or add it to your system PATH.
#
# Input:
#   - rawfile.vcf : raw VCF file from GBS pipeline (e.g. TASSEL or GATK)
#
# Output:
#   - raw.ped / raw.map          : PLINK text format
#   - raw_bed.bed/bim/fam        : PLINK binary format
#   - raw_maf_geno.bed/bim/fam   : QC-filtered binary files
#   - raw_qc_out_vcf.vcf         : QC-filtered VCF for imputation
#
# NOTE — Heterozygosity filtering (TASSEL):
#   Filtering based on 20% heterozygosity per SNP was performed in TASSEL
#   prior to running this script. TASSEL's Filter Genotype Table Sites plugin
#   was used with the following setting:
#     - Maximum proportion of heterozygous sites: 0.20
#   The heterozygosity-filtered VCF was used as the input (rawfile.vcf) here.
#   TASSEL GUI: https://www.maizegenetics.net/tassel
#
# Author:  Rainy
# Contact: rainy122001@gmail.com
# DOI:     https://doi.org/10.1007/s13353-026-01078-3
# License: MIT
# =============================================================================


# =============================================================================
# 0. USER CONFIGURATION — edit this section before running
# =============================================================================

# Path to PLINK executable
# Examples:
#   PLINK_EXE <- "plink"          # Linux/Mac (if added to PATH)
#   PLINK_EXE <- "./plink.exe"    # Windows (in working directory)
PLINK_EXE <- "plink.exe"

# Input VCF file (heterozygosity-filtered in TASSEL)
INPUT_VCF <- "rawfile.vcf"

# Output prefix for each step
OUT_RAW     <- "raw"              # PLINK text format
OUT_BED     <- "raw_bed"          # PLINK binary format
OUT_QC      <- "raw_maf_geno"     # After MAF + missingness filtering
OUT_QC_VCF  <- "raw_qc_out_vcf"  # Final filtered VCF

# QC thresholds
MAF_THRESHOLD  <- 0.05   # Minimum minor allele frequency (5%)
GENO_THRESHOLD <- 0.20   # Maximum missingness per SNP (20%)

# Set working directory — update to your local project folder
setwd("path/to/your/project")


# =============================================================================
# 1. CONVERT VCF TO PLINK TEXT FORMAT
# =============================================================================
# --vcf          : input VCF file
# --recode       : output as PLINK text format (.ped + .map)
# --double-id    : use sample ID for both family and individual ID fields
# --nonfounders  : include all individuals regardless of founder status
# --allow-extra-chr : allow non-standard chromosome names (e.g. maize chr names)
# --allow-no-sex : suppress error when sex information is missing

cat("Step 1: Converting VCF to PLINK text format...\n")

system(paste(
  PLINK_EXE,
  "--vcf", INPUT_VCF,
  "--recode",
  "--double-id",
  "--nonfounders",
  "--allow-extra-chr",
  "--allow-no-sex",
  "--out", OUT_RAW
))

cat("Done. Output:", OUT_RAW, ".ped / .map\n\n")


# =============================================================================
# 2. CONVERT TO PLINK BINARY FORMAT
# =============================================================================
# Binary format (.bed/.bim/.fam) is required for QC filtering steps.
# --file     : input PLINK text format prefix
# --make-bed : output binary format

cat("Step 2: Converting to PLINK binary format...\n")

system(paste(
  PLINK_EXE,
  "--file", OUT_RAW,
  "--make-bed",
  "--double-id",
  "--nonfounders",
  "--allow-extra-chr",
  "--allow-no-sex",
  "--out", OUT_BED
))

cat("Done. Output:", OUT_BED, ".bed/.bim/.fam\n\n")


# =============================================================================
# 3. APPLY QC FILTERS
# =============================================================================
# Filters applied simultaneously:
#   --maf  : removes SNPs with MAF below threshold (rare variants)
#   --geno : removes SNPs with missingness above threshold
#
# NOTE: Individual-level missingness (--mind) was not applied here as
# the panel was curated prior to genotyping. Add --mind 0.20 if needed
# for your dataset.
#
# To customise thresholds, adjust MAF_THRESHOLD and GENO_THRESHOLD in
# Section 0.

cat("Step 3: Applying QC filters (MAF =", MAF_THRESHOLD,
    ", missingness =", GENO_THRESHOLD, ")...\n")

system(paste(
  PLINK_EXE,
  "--bfile", OUT_BED,
  "--allow-extra-chr",
  "--geno", GENO_THRESHOLD,
  "--maf",  MAF_THRESHOLD,
  "--make-bed",
  "--out", OUT_QC
))

cat("Done. Output:", OUT_QC, ".bed/.bim/.fam\n\n")


# =============================================================================
# 4. EXPORT FILTERED GENOTYPES AS VCF
# =============================================================================
# Converts QC-filtered binary files back to VCF format for Beagle imputation.
# --recode vcf-iid : output VCF with individual IDs (not family IDs)

cat("Step 4: Exporting QC-filtered VCF...\n")

system(paste(
  PLINK_EXE,
  "--bfile", OUT_QC,
  "--recode vcf-iid",
  "--allow-extra-chr",
  "--out", OUT_QC_VCF
))

cat("Done. Output:", OUT_QC_VCF, ".vcf\n\n")
cat("=== QC complete. Proceed to 02_imputation.sh ===\n")
cat("Input for imputation:", paste0(OUT_QC_VCF, ".vcf"), "\n")
