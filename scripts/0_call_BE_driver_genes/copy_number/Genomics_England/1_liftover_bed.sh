#!/bin/bash

# Directory containing BED files
BED_DIR="/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/published_version/data/GEL/cna_canvas_pass"

# Directory to output lifted BED files
OUTPUT_DIR="/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/published_version/data/GEL/cna_canvas_pass_hg19"

# Chain file path
CHAIN_FILE="/mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/scripts/picard/hg38ToHg19.over.chain"

# Ensure output directories exist
#mkdir -p "$OUTPUT_DIR" "$UNMAPPED_DIR"

# Loop through each BED file in the directory
for bed_file in "$BED_DIR"/*.bed; do
    echo "Processing $bedfile"
    # Construct output and unmapped file names based on input file
    output_file="$OUTPUT_DIR/$(basename "${bed_file%.bed}")_hg19.bed"
    unmapped_file="$OUTPUT_DIR/$(basename "${bed_file%.bed}")_unmapped.bed"

    # Perform liftover
    /scratchc/stlab-icgc/users/wu04/project/bo_gene_list/scripts/liftover/liftOver "$bed_file" "$CHAIN_FILE" "$output_file" "$unmapped_file"
done

#/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/scripts/liftover/liftOver /scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/gel/cna_canvas_pass_hg19/gel_cna_pass.bed /scratchb/stlab-icgc/users/wu04/project/bo_gene_list/scripts/picard/hg38ToHg19.over.chain /scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/gel/cna_canvas_pass_hg19/gel_cna_pass_hg19.bed /scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/gel/cna_canvas_pass_hg19/gel_cna_pass_unmapped.bed
