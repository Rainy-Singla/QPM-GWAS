#!/bin/bash
# =============================================================================
# Genotype Imputation using Beagle v5.5
# =============================================================================
#
# Description:
#   This script performs imputation of missing genotype data in the
#   QC-filtered VCF file using Beagle v5.5. Beagle uses a hidden Markov
#   model to phase genotypes and impute missing values based on observed
#   haplotype patterns in the dataset.
#
# Prerequisites:
#   - Java (>= 8) installed and accessible from the command line
#   - Beagle JAR file downloaded from:
#     https://faculty.washington.edu/browning/beagle/beagle.html
#
# Input:
#   - QC-filtered VCF from 01_QC.R (raw_qc_out_vcf.vcf)
#
# Output:
#   - Imputed VCF file (imputed_output.vcf.gz)
#   - Unzipped VCF for downstream GWAS (imputed_output.vcf)
#
# Memory:
#   -Xmx8g allocates 8GB of RAM to Java. Increase (e.g. -Xmx16g) if your
#   dataset is large or if the process runs out of memory. Decrease if your
#   system has limited RAM.
#
# Author:  Rainy
# Contact: rainy122001@gmail.com
# DOI:     https://doi.org/10.1007/s13353-026-01078-3
# License: MIT
# =============================================================================


# =============================================================================
# USER CONFIGURATION — edit these variables before running
# =============================================================================

# Path to Beagle JAR file — update to match your downloaded version
BEAGLE_JAR="beagle.27Feb25.75.jar"

# Input VCF — output from 01_QC.R
INPUT_VCF="raw_qc_out_vcf.vcf"

# Output prefix (Beagle appends .vcf.gz automatically)
OUTPUT_PREFIX="imputed_output"

# Memory allocation for Java (in GB)
# Increase if processing large datasets
MEMORY="8g"


# =============================================================================
# IMPUTATION
# =============================================================================
# Parameters:
#   gt  = input VCF file with observed (possibly missing) genotypes
#   out = output file prefix
#
# Optional parameters not used here but available for tuning:
#   nthreads = number of CPU threads (default: all available)
#   window   = length of sliding window in cM (default: 40.0)
#   overlap  = overlap between windows in cM (default: 2.0)
#   iterations = number of phasing iterations (default: 5)
#
# For reference panel-based imputation, add:
#   ref = reference_panel.vcf.gz

echo "Starting Beagle imputation..."
echo "Input:  ${INPUT_VCF}"
echo "Output: ${OUTPUT_PREFIX}.vcf.gz"
echo "Memory: ${MEMORY}"

java -Xmx${MEMORY} -jar ${BEAGLE_JAR} \
    gt=${INPUT_VCF} \
    out=${OUTPUT_PREFIX}

# Check if imputation succeeded
if [ $? -eq 0 ]; then
    echo "Imputation complete. Output: ${OUTPUT_PREFIX}.vcf.gz"
else
    echo "ERROR: Beagle imputation failed. Check Java version and memory allocation."
    exit 1
fi


# =============================================================================
# DECOMPRESS OUTPUT
# =============================================================================
# Beagle outputs a gzipped VCF. Decompress for use in GWAS tools that
# require uncompressed VCF, or keep as .vcf.gz if your tools support it.

echo "Decompressing output..."
gunzip -k ${OUTPUT_PREFIX}.vcf.gz

if [ $? -eq 0 ]; then
    echo "Decompressed: ${OUTPUT_PREFIX}.vcf"
    echo "=== Imputation complete. Proceed to 03_GWAS.R ==="
else
    echo "ERROR: Decompression failed."
    exit 1
fi
