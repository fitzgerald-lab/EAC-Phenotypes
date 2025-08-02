#!/bin/bash

VCF_DIR="/path/to/input/vcf_directory"

# Define output directory for lifted VCFs
OUTPUT_DIR="/path/to/output/vcf_directory"

# Location of the Picard jar
PICARD_JAR="/path/to/picard.jar"

# Chain file for liftover
CHAIN_FILE="/path/to/hg38ToHg19.over.chain"

# Reference sequence fasta
REFERENCE_SEQUENCE="/path/to/hg19.fa"


# Loop through each VCF file in the directory
for VCF in "$VCF_DIR"/*.vcf; do
    OUTPUT_VCF="$OUTPUT_DIR/$(basename "$VCF" .vcf)_lifted_hg19.vcf"
    REJECTED_VCF="$OUTPUT_DIR/$(basename "$VCF" .vcf)_rejected.vcf"

    # Run Picard LiftoverVcf
    java -jar "$PICARD_JAR" LiftoverVcf \
    I="$VCF" \
    O="$OUTPUT_VCF" \
    CHAIN="$CHAIN_FILE" \
    REJECT="$REJECTED_VCF" \
    R="$REFERENCE_SEQUENCE"

    echo "Processed: $VCF"
done
# End of script