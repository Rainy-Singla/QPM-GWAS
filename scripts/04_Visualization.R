# =============================================================================
# GWAS Visualisation — Manhattan Plot and Venn Diagram
# =============================================================================
#
# Description:
#   This script produces publication-ready visualisations from GWAS results:
#
#     1. Manhattan plot — genome-wide p-value distribution with significant
#        SNPs highlighted and labelled (CMplot)
#     2. Venn diagram — overlap of significant SNPs across GWAS models (ggvenn)
#
# Input:
#   - GAPIT GWAS results CSV (all SNPs, one model)
#   - Excel file of selected/stable SNPs to highlight on Manhattan plot
#   - Excel file of significant SNPs per model for Venn diagram
#     (two columns: SNP, Model)
#
# Output:
#   - figures/Manhattan_<trait>_<model>.jpg   — Manhattan plot
#   - figures/Venn_<trait>.png                — Venn diagram
#
# Author:  Rainy
# Contact: rainy122001@gmail.com
# DOI:     https://doi.org/10.1007/s13353-026-01078-3
# License: MIT
# =============================================================================


# =============================================================================
# 0. USER CONFIGURATION — edit this section before running
# =============================================================================

# --- 0.1 Manhattan plot inputs -----------------------------------------------

# GAPIT results CSV for the model and trait to plot
# File format: SNP, Chr, Pos, P.value (standard GAPIT output)
GWAS_CSV <- "results/GAPIT_KNPC/GAPIT.Association.GWAS_Results.FarmCPU.KNPC.csv"

# Excel file containing SNPs to highlight and label on the Manhattan plot
# Expected column: "SNP" (SNP IDs matching those in GWAS_CSV)
# Set to NULL to skip highlighting
HIGHLIGHT_FILE  <- "data/selected_SNP.xlsx"
HIGHLIGHT_SHEET <- 7   # sheet number in the Excel file (ignored if NULL)

# Trait and model name — used for output file naming
TRAIT_NAME <- "KNPC"
MODEL_NAME <- "FarmCPU"

# Significance threshold on the Manhattan plot (as p-value, not -log10)
# 10^-3.5 ~ 3.2e-4; adjust based on your multiple testing approach
SIGNIFICANCE_THRESHOLD <- 10^-3.5

# --- 0.2 Venn diagram inputs -------------------------------------------------

# Excel file with two columns: SNP and Model
# Each row = one significant SNP from one model
# Models column should match model names exactly (e.g. "FarmCPU", "MLM")
VENN_FILE <- "results/significant_SNPs_ALL_models_KNPC.xlsx"

# Models to include in the Venn diagram
# Must match values in the "Model" column of VENN_FILE
VENN_MODELS <- c("FarmCPU", "MLM", "BLINK", "CMLM")

# Colours for Venn diagram circles (one per model)
VENN_COLOURS <- c("skyblue", "lightgreen", "lightpink", "orchid")

# Set working directory — update to your local project folder
# setwd("path/to/your/project")


# =============================================================================
# 1. SETUP
# =============================================================================

required_packages <- c("CMplot", "ggvenn", "readxl", "dplyr", "ggplot2")

new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages) > 0) {
  message("Installing missing packages: ", paste(new_packages, collapse = ", "))
  install.packages(new_packages)
}

suppressPackageStartupMessages({
  library(CMplot)
  library(ggvenn)
  library(readxl)
  library(dplyr)
  library(ggplot2)
})

dir.create("figures", showWarnings = FALSE)

cat("=== GWAS Visualisation ===\n")
cat("Trait:", TRAIT_NAME, "| Model:", MODEL_NAME, "\n\n")


# =============================================================================
# 2. MANHATTAN PLOT
# =============================================================================
# Plots genome-wide -log10(p) values per chromosome.
# Significant SNPs are highlighted and labelled.
#
# CMplot parameters used:
#   plot.type = "m"  : Manhattan plot
#   LOG10     = TRUE : converts p-values to -log10 scale
#   highlight        : SNP IDs to colour differently
#   highlight.text   : labels shown above highlighted SNPs
#
# NOTE: P-values of exactly 0 are replaced with 1e-300 to avoid -log10(0) = Inf.
# P-values of exactly 1 are replaced with 0.999999999 for display purposes.

cat("--- Section 2: Manhattan Plot ---\n")

# Load all GWAS results
gwas_results <- read.csv(GWAS_CSV, header = TRUE)

# Keep only required columns and rename for CMplot
gwas_results <- gwas_results[, c("SNP", "Chr", "Pos", "P.value")]
colnames(gwas_results) <- c("SNP", "CHR", "BP", "P")
gwas_results$P <- as.numeric(gwas_results$P)

# Replace boundary p-values
gwas_results$P[gwas_results$P == 0] <- 1e-300
gwas_results$P[gwas_results$P >= 1] <- 0.999999999

cat("  Total SNPs:", nrow(gwas_results), "\n")

# Load SNPs to highlight
if (!is.null(HIGHLIGHT_FILE) && file.exists(HIGHLIGHT_FILE)) {
  selected_snps <- read_excel(HIGHLIGHT_FILE, sheet = HIGHLIGHT_SHEET)
  label_SNP     <- selected_snps$SNP
  cat("  SNPs to highlight:", length(label_SNP), "\n")
} else {
  label_SNP <- NULL
  cat("  No highlight file provided — plotting without highlights\n")
}

# Generate Manhattan plot
CMplot(
  gwas_results,
  plot.type           = "m",
  LOG10               = TRUE,
  threshold           = SIGNIFICANCE_THRESHOLD,
  threshold.lty       = 2,
  threshold.col       = "black",
  threshold.lwd       = 2,
  highlight           = label_SNP,
  highlight.col       = "black",
  highlight.text      = label_SNP,
  highlight.text.col  = "black",
  highlight.text.cex  = 1.5,
  chr.labels.angle    = 0,
  col                 = c("red", "orange", "#66C266",
                           "steelblue", "#A366CC"),
  file                = "jpg",
  file.name           = paste0("Manhattan_", TRAIT_NAME, "_", MODEL_NAME),
  file.output         = TRUE,
  output.file         = paste0("figures/Manhattan_", TRAIT_NAME, "_", MODEL_NAME),
  dpi                 = 300
)

cat("  Saved: figures/Manhattan_", TRAIT_NAME, "_", MODEL_NAME, ".jpg\n\n",
    sep = "")


# =============================================================================
# 3. VENN DIAGRAM
# =============================================================================
# Shows the overlap of significant SNPs identified across GWAS models.
# SNPs shared between models are more likely to be true positives.
#
# The diagram uses ggvenn which supports up to 4 sets.
# If you have more than 4 models, consider using UpSetR instead
# (see commented alternative below).

cat("--- Section 3: Venn Diagram ---\n")

# Load significant SNPs file
if (!file.exists(VENN_FILE)) {
  stop("Venn diagram input file not found: ", VENN_FILE,
       "\nRun 03_GWAS.R first to generate significant SNP files.")
}

snp_data <- read_excel(VENN_FILE)

# Validate expected columns
if (!all(c("SNP", "Model") %in% colnames(snp_data))) {
  stop("VENN_FILE must contain columns named 'SNP' and 'Model'.")
}

# Filter to specified models and remove duplicates
snp_data <- snp_data %>%
  filter(Model %in% VENN_MODELS) %>%
  distinct(SNP, Model)

# Build named list of SNP sets per model
snp_sets <- split(snp_data$SNP, snp_data$Model)
snp_sets <- lapply(snp_sets, unique)

cat("  SNPs per model:\n")
for (m in names(snp_sets)) {
  cat("   ", m, ":", length(snp_sets[[m]]), "\n")
}

# Validate colour length
if (length(VENN_COLOURS) < length(snp_sets)) {
  warning("Fewer colours than models — recycling colours.")
  VENN_COLOURS <- rep(VENN_COLOURS, length.out = length(snp_sets))
}

# Generate Venn diagram
p_venn <- ggvenn(
  snp_sets,
  fill_color   = VENN_COLOURS[seq_along(snp_sets)],
  stroke_size  = 0.8,
  set_name_size = 5,
  text_size    = 4
)

venn_out <- paste0("figures/Venn_", TRAIT_NAME, ".png")
ggsave(venn_out, plot = p_venn, width = 8, height = 7, dpi = 300)
cat("  Saved:", venn_out, "\n\n")


# =============================================================================
# ALTERNATIVE: UpSet plot (for > 4 models or complex overlaps)
# =============================================================================
# Uncomment and run if you prefer an UpSet plot over a Venn diagram.
# UpSet plots handle more than 4 sets and show intersection sizes clearly.
#
# if (!requireNamespace("UpSetR", quietly = TRUE))
#   install.packages("UpSetR")
# library(UpSetR)
#
# png(paste0("figures/UpSet_", TRAIT_NAME, ".png"),
#     width = 3000, height = 2000, res = 300)
# upset(fromList(snp_sets),
#       nsets        = length(snp_sets),
#       order.by     = "freq",
#       sets.bar.color = "steelblue")
# dev.off()


# =============================================================================
# END
# =============================================================================

cat("=== Visualisation complete ===\n")
cat("Figures saved to: figures/\n")
sessionInfo()
