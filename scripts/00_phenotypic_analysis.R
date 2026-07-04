# =============================================================================
# Phenotypic Analysis Pipeline for Resolvable Incomplete Block (Lattice) Designs
# =============================================================================
#
# Description:
#   A generalised, reusable pipeline for the statistical analysis of
#   multi-trait phenotypic data from plant breeding trials conducted in a
#   resolvable incomplete block (lattice) design. The script covers:
#
#     1.  Data import         — Excel (.xlsx) or CSV
#     2.  Descriptive stats   — mean, median, SD, CV%, skewness, kurtosis
#     3.  Boxplots            — compact scattered boxplot per trait
#     4.  Frequency plots     — histogram + density per trait
#     5.  ANOVA               — lattice model with outlier detection
#     6.  Critical Difference — CD at 5% significance
#     7.  Genetic variability — GV, PV, h², GA, GA%
#     8.  Normality testing   — Shapiro-Wilk + Box-Cox transformation
#     9.  BLUE estimation     — via mixed model (lme4 + emmeans)
#     10. Correlation         — Pearson matrix, corrplot, and ggpairs panel
#
# Experimental design supported:
#   Resolvable incomplete block (lattice) design
#   Model: trait ~ Genotype + Replication + Replication:Block
#
# How to use this script:
#   1. Fill in SECTION 0 — USER CONFIGURATION below
#   2. Ensure your input file columns match the names you define there
#   3. Run the script section by section, or source() it all at once
#   4. All outputs are written to results/ and figures/ automatically
#
# Input file format:
#   One row per plot. Required columns:
#     - Genotype column    (name defined in COL_GENOTYPE)
#     - Replication column (name defined in COL_REPLICATION)
#     - Block column       (name defined in COL_BLOCK)
#     - One or more numeric trait columns
#
# Outputs:
#   results/Descriptive_Stats.xlsx   — descriptive statistics and CV table
#   results/ANOVA_CD.xlsx            — ANOVA table and CD values
#   results/Genetic_Variability.xlsx — GV, PV, h², GA, GA%
#   results/BLUEs_Normality.xlsx     — BLUEs and normality check
#   results/Correlation.xlsx         — correlation matrix and p-values
#   figures/                         — all plots as high-resolution PNGs/TIFFs
#
# Author:  Rainy
# Contact: rainy122001@gmail.com
# DOI:     https://doi.org/10.1007/s13353-026-01078-3
# License: MIT
# =============================================================================


# =============================================================================
# 0. USER CONFIGURATION — edit only this section before running
# =============================================================================

# --- 0.1 Input file ----------------------------------------------------------

# Full or relative path to your data file (Excel or CSV)
# Examples:
#   INPUT_FILE <- "data/trial_data.xlsx"
#   INPUT_FILE <- "data/trial_data.csv"
INPUT_FILE <- "data/Phenotypic_data.xlsx"

# For Excel only: which sheet number contains the raw data? (ignored for CSV)
SHEET <- 2

# --- 0.2 Column names for design factors ------------------------------------
# Set these to match your actual column names exactly (case-sensitive)

COL_GENOTYPE    <- "TREATMENT"    # Column identifying genotypes/inbred lines
COL_REPLICATION <- "REPLICATION"  # Column identifying replications
COL_BLOCK       <- "BLOCK"        # Column identifying incomplete blocks

# --- 0.3 Trait columns -------------------------------------------------------
# Option A — Manual: list your trait column names explicitly
#   Example: TRAITS_MANUAL <- c("KL", "KW", "KT", "CD")
#   Leave as NULL to use auto-detection (Option B)
#
# Option B — Auto-detect: identifies all numeric columns that are not
#   design factor columns. Detected traits are printed at runtime for review.

TRAITS_MANUAL <- NULL   # e.g. c("KL", "KW", "KT", "CD", "CL", "KNPR", "RN", "KNPC")
                        # Set to NULL for auto-detection

# --- 0.4 Experimental design parameters -------------------------------------

N_REPS <- 2   # Number of replications per environment

# --- 0.5 Output file paths --------------------------------------------------
# Rename if needed; all files are saved inside results/ and figures/ folders

OUT_DESC  <- "results/Descriptive_Stats.xlsx"
OUT_ANOVA <- "results/ANOVA_CD.xlsx"
OUT_GV    <- "results/Genetic_Variability.xlsx"
OUT_BLUE  <- "results/BLUEs_Normality.xlsx"
OUT_COR   <- "results/Correlation.xlsx"


# =============================================================================
# 1. SETUP — install missing packages and load libraries
# =============================================================================

required_packages <- c(
  "dplyr", "ggplot2", "emmeans", "openxlsx", "readxl",
  "tidyr", "moments", "gridExtra", "purrr", "MASS",
  "lme4", "psych", "tibble", "Hmisc", "corrplot",
  "GGally", "scales", "cowplot", "ggcorrplot"
)

new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages) > 0) {
  message("Installing missing packages: ", paste(new_packages, collapse = ", "))
  install.packages(new_packages)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(emmeans)
  library(openxlsx)
  library(readxl)
  library(tidyr)
  library(moments)
  library(gridExtra)
  library(purrr)
  library(MASS)
  library(lme4)
  library(psych)
  library(tibble)
  library(Hmisc)
  library(corrplot)
  library(GGally)
  library(scales)
  library(cowplot)
})

# Create output directories if they do not already exist
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cat("=== Phenotypic Analysis Pipeline ===\n")
cat("R version:", paste0(R.version$major, ".", R.version$minor), "\n")


# =============================================================================
# 2. DATA IMPORT AND VALIDATION
# =============================================================================

file_ext <- tolower(tools::file_ext(INPUT_FILE))

if (file_ext %in% c("xlsx", "xls")) {
  data <- read_excel(INPUT_FILE, sheet = SHEET)
  cat("Imported from Excel:", INPUT_FILE, "(sheet", SHEET, ")\n")
} else if (file_ext == "csv") {
  data <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
  cat("Imported from CSV:", INPUT_FILE, "\n")
} else {
  stop("Unsupported file type '", file_ext,
       "'. Please provide an .xlsx or .csv file.")
}

# Validate that required design columns exist
required_cols <- c(COL_GENOTYPE, COL_REPLICATION, COL_BLOCK)
missing_cols  <- required_cols[!(required_cols %in% colnames(data))]
if (length(missing_cols) > 0) {
  stop(
    "Required column(s) not found in data:\n  ",
    paste(missing_cols, collapse = ", "),
    "\nCheck COL_GENOTYPE, COL_REPLICATION, and COL_BLOCK in Section 0."
  )
}

# Standardise design factor names internally
data$Genotype    <- as.factor(data[[COL_GENOTYPE]])
data$Replication <- as.factor(data[[COL_REPLICATION]])
data$Block       <- as.factor(data[[COL_BLOCK]])

# --- Trait detection --------------------------------------------------------
if (!is.null(TRAITS_MANUAL)) {
  missing_traits <- TRAITS_MANUAL[!(TRAITS_MANUAL %in% colnames(data))]
  if (length(missing_traits) > 0)
    stop("Manually specified trait(s) not found:\n  ",
         paste(missing_traits, collapse = ", "))
  traits <- TRAITS_MANUAL
  cat("Traits (manual):", paste(traits, collapse = ", "), "\n")
} else {
  factor_cols <- c(COL_GENOTYPE, COL_REPLICATION, COL_BLOCK,
                   "Genotype", "Replication", "Block")
  traits <- colnames(data)[
    sapply(data, is.numeric) & !(colnames(data) %in% factor_cols)
  ]
  if (length(traits) == 0)
    stop("No numeric trait columns detected. ",
         "Set TRAITS_MANUAL in Section 0.")
  cat("Traits (auto-detected):", paste(traits, collapse = ", "), "\n")
  cat("If incorrect, set TRAITS_MANUAL manually in Section 0.\n")
}

cat("\nDataset overview:\n")
cat("  Rows:", nrow(data), "\n")
cat("  Genotypes:", nlevels(data$Genotype), "\n")
cat("  Replications:", nlevels(data$Replication), "\n")
cat("  Traits:", length(traits), "\n\n")


# =============================================================================
# 3. DESCRIPTIVE STATISTICS
# =============================================================================
# Per-trait summary: mean, median, SD, min, max, skewness, and kurtosis.
# Skewness and kurtosis significance tested via Z-score:
#   SE_skewness = sqrt(6/n),  SE_kurtosis = sqrt(24/n)

cat("--- Section 3: Descriptive Statistics ---\n")

desc_stats <- data %>%
  summarise(across(all_of(traits), list(
    Mean     = ~ mean(.x,     na.rm = TRUE),
    Median   = ~ median(.x,   na.rm = TRUE),
    SD       = ~ sd(.x,       na.rm = TRUE),
    Min      = ~ min(.x,      na.rm = TRUE),
    Max      = ~ max(.x,      na.rm = TRUE),
    Skewness = ~ skewness(.x, na.rm = TRUE),
    Kurtosis = ~ kurtosis(.x, na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(),
               names_to  = c("Trait", "Statistic"),
               names_sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = "Statistic", values_from = "value")

# Z-score significance tests for skewness and kurtosis
skew_kurt_p <- map_dfr(traits, function(trait) {
  x      <- na.omit(data[[trait]])
  n      <- length(x)
  skew   <- skewness(x)
  kurt   <- kurtosis(x)
  z_skew <- skew / sqrt(6 / n)
  z_kurt <- (kurt - 3) / sqrt(24 / n)
  data.frame(
    Trait  = trait,
    p_skew = round(2 * (1 - pnorm(abs(z_skew))), 4),
    p_kurt = round(2 * (1 - pnorm(abs(z_kurt))), 4)
  )
})

desc_stats <- left_join(desc_stats, skew_kurt_p, by = "Trait")
print(desc_stats)

# --- Coefficient of Variation -----------------------------------------------
cv_table <- map_dfr(traits, function(trait) {
  m <- mean(data[[trait]], na.rm = TRUE)
  s <- sd(data[[trait]],   na.rm = TRUE)
  data.frame(
    Trait             = trait,
    Mean              = round(m, 3),
    SD                = round(s, 3),
    CV_percent        = round((s / m) * 100, 2),
    Variability_level = case_when(
      (s / m) * 100 < 10 ~ "Low",
      (s / m) * 100 < 20 ~ "Moderate",
      TRUE               ~ "High"
    )
  )
})
print(cv_table)

# Save descriptive stats
wb_desc <- createWorkbook()
addWorksheet(wb_desc, "Descriptive_Stats")
writeData(wb_desc, "Descriptive_Stats", desc_stats)
addWorksheet(wb_desc, "CV_Table")
writeData(wb_desc, "CV_Table", cv_table)
saveWorkbook(wb_desc, OUT_DESC, overwrite = TRUE)
cat("Saved:", OUT_DESC, "\n")


# =============================================================================
# 4. BOXPLOTS — compact scattered boxplot per trait
# =============================================================================
# Displays both the distributional summary (box) and individual data points
# (jittered dots) simultaneously. Outliers shown as jittered points rather
# than as separate symbols, avoiding double-plotting.

cat("--- Section 4: Boxplots ---\n")

walk(traits, function(trait) {
  p <- ggplot(data, aes(x = "", y = .data[[trait]])) +
    geom_boxplot(
      width        = 0.3,
      fill         = "#7FB3D5",
      color        = "black",
      outlier.shape = NA          # hide default outlier points; shown via jitter
    ) +
    geom_jitter(
      width = 0.1,
      alpha = 0.6,
      color = "#2E86C1",
      size  = 1.8
    ) +
    labs(
      title = paste("Distribution of", trait),
      y     = paste(trait, "value"),
      x     = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )

  ggsave(
    filename = paste0("figures/Boxplot_", trait, ".png"),
    plot     = p,
    width    = 4, height = 5, dpi = 300
  )
  cat("  Saved: figures/Boxplot_", trait, ".png\n", sep = "")
})


# =============================================================================
# 5. FREQUENCY DISTRIBUTIONS — histogram + density per trait
# =============================================================================
# Red solid line = mean; green dashed line = median.
# Helps assess normality and distribution shape before formal testing.

cat("--- Section 5: Frequency Distributions ---\n")

walk(traits, function(trait) {
  x          <- na.omit(data[[trait]])
  mean_val   <- mean(x)
  median_val <- median(x)

  p <- ggplot(data.frame(x = x), aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins  = 15, fill = "grey80", color = "black") +
    geom_density(color = "royalblue", linewidth = 0.9) +
    geom_vline(xintercept = mean_val,   color = "red",
               linetype = "solid",  linewidth = 0.9) +
    geom_vline(xintercept = median_val, color = "darkgreen",
               linetype = "dashed", linewidth = 0.9) +
    annotate("text", x = mean_val,   y = Inf,
             label = "Mean",   color = "red",
             vjust = 2, hjust = -0.15, size = 3.5) +
    annotate("text", x = median_val, y = Inf,
             label = "Median", color = "darkgreen",
             vjust = 2, hjust = -0.15, size = 3.5) +
    labs(title = paste("Frequency Distribution —", trait),
         x = trait, y = "Density") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(
    filename = paste0("figures/FreqDist_", trait, ".png"),
    plot     = p,
    width = 5, height = 4, dpi = 300
  )
  cat("  Saved: figures/FreqDist_", trait, ".png\n", sep = "")
})


# =============================================================================
# 6. ANOVA — resolvable incomplete block (lattice) design
# =============================================================================
# Model: trait ~ Genotype + Replication + Replication:Block
#
# Sources of variation (SV):
#   Genotype            — genetic differences among lines
#   Replication         — systematic differences between replications
#   Block:Replication   — incomplete block effect within replication
#   Residuals           — unexplained experimental error
#
# Outlier detection: plots with |studentized residual| > 3 are removed
# before the final model is fitted.

cat("--- Section 6: ANOVA ---\n")

anova_results <- map_dfr(seq_along(traits), function(i) {
  trait <- traits[i]
  cat("  Trait", i, "of", length(traits), ":", trait, "\n")

  form <- as.formula(
    paste(trait, "~ Genotype + Replication + Replication:Block")
  )

  # Outlier removal
  fit_init    <- lm(form, data = data)
  outlier_idx <- which(abs(rstudent(fit_init)) > 3)
  data_clean  <- if (length(outlier_idx) == 0) data else data[-outlier_idx, ]
  if (length(outlier_idx) > 0)
    cat("    Outliers removed:", length(outlier_idx), "\n")

  fit_final <- lm(form, data = data_clean)
  anv       <- as.data.frame(anova(fit_final))
  anv$Trait <- trait
  anv       <- add_column(anv,
    SV = c("Genotype", "Replication", "Block:Replication", "Residuals"),
    .before = "Df"
  )
  rownames(anv) <- NULL
  anv
})

wb_anova <- createWorkbook()
addWorksheet(wb_anova, "ANOVA_Table")
writeData(wb_anova, "ANOVA_Table", anova_results)
cat("ANOVA complete.\n")


# =============================================================================
# 7. CRITICAL DIFFERENCE (CD) AT 5%
# =============================================================================
# CD = t(alpha/2, df_error) * sqrt(2 * MSE / r)
# Minimum difference between two genotype means to be declared significant.

cat("--- Section 7: Critical Difference ---\n")

cd_table <- map_dfr(traits, function(trait) {
  model  <- lm(
    as.formula(paste(trait, "~ Genotype + Replication + Replication:Block")),
    data = data
  )
  anv    <- anova(model)
  mse    <- tail(anv$`Mean Sq`, 1)
  df_res <- tail(anv$Df,        1)
  t_val  <- qt(0.975, df = df_res)
  cd     <- t_val * sqrt(2 * mse / N_REPS)
  data.frame(
    Trait   = trait,
    MSE     = round(mse, 4),
    CD_5pct = round(cd,  4)
  )
})

print(cd_table)

addWorksheet(wb_anova, "CD_Table")
writeData(wb_anova, "CD_Table", cd_table)
saveWorkbook(wb_anova, OUT_ANOVA, overwrite = TRUE)
cat("Saved:", OUT_ANOVA, "\n")


# =============================================================================
# 8. GENETIC VARIABILITY PARAMETERS
# =============================================================================
# Variance components estimated from ANOVA mean squares:
#   GV  = Genotypic Variance        = (MSG - MSE) / r    [floored at 0]
#   PV  = Phenotypic Variance       = GV + MSE
#   h²  = Broad-sense Heritability  = GV / PV
#   GA  = Genetic Advance           = 2.06 * sqrt(PV) * h²  [5% selection]
#   GA% = GA as % of grand mean
#
# Category thresholds (Burton & DeVane, 1953):
#   h²:  High >= 0.60 | Moderate 0.30-0.59 | Low < 0.30
#   GA%: High >= 20%  | Moderate 10-19%    | Low < 10%

cat("--- Section 8: Genetic Variability ---\n")

gv_stats <- map_dfr(traits, function(trait) {
  model   <- aov(
    as.formula(paste(trait, "~ Genotype + Replication + Replication:Block")),
    data = data
  )
  aov_tbl    <- anova(model)
  MSG        <- aov_tbl$`Mean Sq`[1]
  MSE        <- aov_tbl$`Mean Sq`[4]
  GV         <- max((MSG - MSE) / N_REPS, 0)
  PV         <- GV + MSE
  h2         <- ifelse(PV > 0, GV / PV, NA)
  GA         <- 2.06 * sqrt(PV) * h2
  grand_mean <- mean(data[[trait]], na.rm = TRUE)
  GA_pct     <- (GA / grand_mean) * 100

  data.frame(
    Trait               = trait,
    Grand_Mean          = round(grand_mean, 3),
    Genotypic_Variance  = round(GV,     4),
    Phenotypic_Variance = round(PV,     4),
    Heritability_h2     = round(h2,     3),
    Genetic_Advance     = round(GA,     3),
    GA_percent_of_mean  = round(GA_pct, 2),
    h2_category = case_when(
      h2 >= 0.60 ~ "High", h2 >= 0.30 ~ "Moderate", TRUE ~ "Low"),
    GA_category = case_when(
      GA_pct >= 20 ~ "High", GA_pct >= 10 ~ "Moderate", TRUE ~ "Low")
  )
})

print(gv_stats)

# Heritability bar plot
p_h2 <- ggplot(gv_stats,
               aes(x = reorder(Trait, -Heritability_h2),
                   y = Heritability_h2, fill = h2_category)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  scale_fill_manual(values = c("High"     = "#2E8B57",
                                "Moderate" = "#FFA500",
                                "Low"      = "#CD5C5C")) +
  geom_hline(yintercept = c(0.30, 0.60),
             linetype = "dashed", color = "grey40", linewidth = 0.7) +
  labs(title = "Broad-Sense Heritability (h²)",
       y = "h²", x = NULL, fill = "Category") +
  ylim(0, 1) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold", hjust = 0.5))

# Genetic advance bar plot
p_ga <- ggplot(gv_stats,
               aes(x = reorder(Trait, -GA_percent_of_mean),
                   y = GA_percent_of_mean, fill = GA_category)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  scale_fill_manual(values = c("High"     = "#FF8C00",
                                "Moderate" = "#87CEEB",
                                "Low"      = "#D3D3D3")) +
  geom_hline(yintercept = c(10, 20),
             linetype = "dashed", color = "grey40", linewidth = 0.7) +
  labs(title = "Genetic Advance as % of Mean",
       y = "GA (%)", x = NULL, fill = "Category") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold", hjust = 0.5))

png("figures/Genetic_Variability.png", width = 3200, height = 1800, res = 300)
grid.arrange(p_h2, p_ga, ncol = 2)
dev.off()
cat("Saved: figures/Genetic_Variability.png\n")

wb_gv <- createWorkbook()
addWorksheet(wb_gv, "Genetic_Variability")
writeData(wb_gv, "Genetic_Variability", gv_stats)
saveWorkbook(wb_gv, OUT_GV, overwrite = TRUE)
cat("Saved:", OUT_GV, "\n")


# =============================================================================
# 9. NORMALITY TESTING, BOX-COX TRANSFORMATION, AND BLUE ESTIMATION
# =============================================================================
# Workflow per trait:
#   Step 1 — Shapiro-Wilk test (suitable for n < 5000)
#   Step 2 — If p < 0.05: apply Box-Cox transformation (optimal lambda)
#             If p >= 0.05: proceed without transformation
#   Step 3 — Missing values imputed with trait-wise mean
#   Step 4 — Mixed model fitted:
#             trait ~ Genotype + Replication + (1|Replication:Block)
#             Genotype and Replication = fixed effects
#             Block:Replication = random effect
#   Step 5 — BLUEs extracted via emmeans()
#
# BLUEs represent genotype means adjusted for replication and block effects
# and are the recommended phenotypic input for GWAS.

cat("--- Section 9: Normality, Transformation, and BLUEs ---\n")

# Impute missing values with trait-wise means
data_imputed <- data
for (trait in traits) {
  na_idx <- is.na(data_imputed[[trait]])
  if (any(na_idx)) {
    data_imputed[[trait]][na_idx] <- mean(data_imputed[[trait]], na.rm = TRUE)
    cat("  Imputed", sum(na_idx), "missing values for:", trait, "\n")
  }
}

norm_results     <- data.frame()
blue             <- data.frame(Genotype = levels(data$Genotype))
transformed_data <- data_imputed

for (trait in traits) {
  cat("  Processing:", trait, "\n")

  sw    <- shapiro.test(data_imputed[[trait]])
  trans <- "None"

  if (sw$p.value < 0.05) {
    bc     <- boxcox(lm(data_imputed[[trait]] ~ 1), plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]
    trans  <- paste0("Box-Cox (lambda = ", round(lambda, 3), ")")

    transformed_data[[trait]] <- if (abs(lambda) < 1e-6) {
      log(data_imputed[[trait]])
    } else {
      (data_imputed[[trait]]^lambda - 1) / lambda
    }
    cat("    Non-normal (p =", round(sw$p.value, 4), ") ->", trans, "\n")
  } else {
    cat("    Normal (p =", round(sw$p.value, 4), ") -> No transformation\n")
  }

  norm_results <- rbind(norm_results, data.frame(
    Trait          = trait,
    Shapiro_W      = round(sw$statistic, 4),
    p_value        = round(sw$p.value,   4),
    Normal         = sw$p.value >= 0.05,
    Transformation = trans
  ))

  # Fit mixed model
  model_mm <- lmer(
    as.formula(paste(trait,
      "~ Genotype + Replication + (1|Replication:Block)")),
    data    = transformed_data,
    control = lmerControl(optimizer = "bobyqa")
  )

  # Extract BLUEs
  em            <- emmeans(model_mm, specs = ~ Genotype)
  blue[[trait]] <- as.numeric(summary(em)$emmean)
}

cat("\nNormality summary:\n")
print(norm_results)

wb_blue <- createWorkbook()
addWorksheet(wb_blue, "BLUEs")
writeData(wb_blue, "BLUEs", blue)
addWorksheet(wb_blue, "Normality_Check")
writeData(wb_blue, "Normality_Check", norm_results)
saveWorkbook(wb_blue, OUT_BLUE, overwrite = TRUE)
cat("Saved:", OUT_BLUE, "\n")


# =============================================================================
# 10. PEARSON CORRELATION ANALYSIS
# =============================================================================
# Computed on BLUEs to reflect genetic correlations among traits,
# free from environmental and block confounding effects.
#
# Three complementary visualisations are produced:
#   A. corrplot   — upper-triangle heatmap, hierarchically clustered,
#                   non-significant pairs (p >= 0.05) shown as blank
#   B. ggpairs    — full pairwise panel: scatter plots (lower), histograms
#                   (diagonal), colour-coded correlation tiles (upper)
#                   with a standalone legend

cat("--- Section 10: Correlation Analysis ---\n")

# Prepare BLUE trait matrix — remove zero-variance columns (safety check)
blue_traits <- blue[, traits, drop = FALSE]
blue_traits <- blue_traits[,
  apply(blue_traits, 2, function(x) sd(x, na.rm = TRUE) > 0),
  drop = FALSE
]
active_traits <- colnames(blue_traits)

# Convert to matrix and remove incomplete rows
traits_matrix <- as.matrix(blue_traits)
traits_matrix <- traits_matrix[
  complete.cases(traits_matrix) &
  apply(traits_matrix, 1, function(x) all(is.finite(x))), ,
  drop = FALSE
]

# Pearson correlation and p-values via Hmisc::rcorr()
cor_results <- rcorr(traits_matrix, type = "pearson")
cat("Correlation matrix:\n")
print(round(cor_results$r, 3))
cat("\nP-values:\n")
print(round(cor_results$P, 4))

# Save correlation outputs
wb_cor <- createWorkbook()
addWorksheet(wb_cor, "Correlation_Matrix")
writeData(wb_cor, "Correlation_Matrix",
  cbind(Trait = rownames(cor_results$r),
        as.data.frame(round(cor_results$r, 4))))
addWorksheet(wb_cor, "P_values")
writeData(wb_cor, "P_values",
  cbind(Trait = rownames(cor_results$P),
        as.data.frame(round(cor_results$P, 4))))
saveWorkbook(wb_cor, OUT_COR, overwrite = TRUE)
cat("Saved:", OUT_COR, "\n")

# --- 10A. corrplot — upper-triangle heatmap ----------------------------------
# Hierarchical clustering reorders traits by similarity.
# Non-significant pairs (p >= 0.05) are left blank.

png("figures/Correlation_corrplot.png", width = 1800, height = 1600, res = 300)
corrplot(
  round(cor_results$r, 2),
  method     = "color",
  type       = "upper",
  order      = "hclust",
  addCoef.col = "black",
  tl.col     = "black",
  tl.srt     = 45,
  p.mat      = cor_results$P,
  sig.level  = 0.05,
  insig      = "blank"
)
dev.off()
cat("Saved: figures/Correlation_corrplot.png\n")

# --- 10B. ggpairs — full pairwise panel with standalone legend ---------------

# Diagonal: histogram of trait values
custom_hist <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_histogram(bins = 20, fill = "#6D9EC1", color = "black") +
    theme_minimal()
}

# Upper panel: colour-coded correlation tile with r value
colored_cor <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  r <- cor(x, y, use = "pairwise.complete.obs")

  df <- data.frame(x = 0.5, y = 0.5, r = ifelse(is.na(r), 0, r))

  ggplot(df, aes(x, y, fill = r)) +
    geom_tile() +
    geom_text(
      aes(label = formatC(r, digits = 2, format = "f")),
      size = 5, color = "black"
    ) +
    scale_fill_gradient2(
      low      = "#B2182B",
      mid      = "white",
      high     = "#2166AC",
      midpoint = 0,
      limits   = c(-1, 1)
    ) +
    theme_void() +
    theme(legend.position = "none")
}

# Build ggpairs panel
p_pairs <- ggpairs(
  as.data.frame(traits_matrix),
  lower = list(continuous = wrap("smooth", alpha = 0.4, color = "#4C4C4C")),
  diag  = list(continuous = custom_hist),
  upper = list(continuous = colored_cor)
)

# Standalone colour legend
legend_plot <- ggplot(
  data.frame(r = seq(-1, 1, length.out = 100)),
  aes(x = r, y = 1, fill = r)
) +
  geom_tile() +
  scale_fill_gradient2(
    low      = "#B2182B",
    mid      = "white",
    high     = "#2166AC",
    midpoint = 0,
    limits   = c(-1, 1),
    name     = "Pearson r"
  ) +
  theme_void() +
  theme(legend.position = "right")

legend_grob <- cowplot::get_legend(legend_plot)

# Combine panel + legend
final_pairs <- cowplot::plot_grid(p_pairs, legend_grob, rel_widths = c(1, 0.08))

ggsave("figures/Correlation_ggpairs.tiff",
       plot   = final_pairs,
       width  = 11, height = 10, dpi = 300)
cat("Saved: figures/Correlation_ggpairs.tiff\n")


# =============================================================================
# END OF PIPELINE
# =============================================================================

cat("\n=== Analysis complete ===\n")
cat("Results: results/\n")
cat("Figures: figures/\n\n")
sessionInfo()
