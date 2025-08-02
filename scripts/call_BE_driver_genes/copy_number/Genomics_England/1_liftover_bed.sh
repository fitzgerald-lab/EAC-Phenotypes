#!/bin/bash

# Directory containing BED files
BED_DIR="./Genomic_analysis/data/0_call_BE_driver_genes/copy_number/Genomics_England"

# Directory to output lifted BED files
OUTPUT_DIR="./Genomic_analysis/data/0_call_BE_driver_genes/copy_number/Genomics_England/lifted_hg19"

# Chain file path
CHAIN_FILE="./Genomic_analysis/scripts/0_call_BE_driver_genes/copy_number/Genomics_England/picard/hg38ToHg19.over.chain"

# Ensure output directories exist
#mkdir -p "$OUTPUT_DIR" "$UNMAPPED_DIR"

# Loop through each BED file in the directory
for bed_file in "$BED_DIR"/*.bed; do
    echo "Processing $bed_file"

    # Remove the first column and header, save to a temporary file
    cleaned_bed="/tmp/cleaned_$(basename "$bed_file")"
    tail -n +2 "$bed_file" | cut -f2- | \
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4":"$5":"$6":"$7":"$8}' > "$cleaned_bed"

    # head "$cleaned_bed"

    # Paths
    lifted_raw="$OUTPUT_DIR/$(basename "${bed_file%.bed}")_hg19_raw.bed"
    lifted_final="$OUTPUT_DIR/$(basename "${bed_file%.bed}")_hg19.bed"
    unmapped_file="$OUTPUT_DIR/$(basename "${bed_file%.bed}")_unmapped.bed"

    # Liftover
    /scratchc/stlab-icgc/users/wu04/project/bo_gene_list/scripts/liftover/liftOver \
        -multiple \
        "$cleaned_bed" "$CHAIN_FILE" "$lifted_raw" "$unmapped_file"

    # Extract original extra columns from name field
    awk 'BEGIN{OFS="\t"} {
        split($4, a, ":");
        $4=a[1]; $5=a[2]; $6=a[3]; $7=a[4]; $8=a[5];
        print $0
    }' "$lifted_raw" > "$lifted_final"

    echo "Saved: $lifted_final"
done
