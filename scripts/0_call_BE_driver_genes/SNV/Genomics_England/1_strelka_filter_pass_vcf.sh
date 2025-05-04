#!/bin/bash

# Input =  ../snv_hg19/*.vcf
# Output =  ../snv_hg19_filter_pass/*_pass.vcf

# Directory containing your VCF files
VCF_DIR="/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/gel/snv_hg19"
 
# Output directory for filtered VCF files
OUTPUT_DIR="/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/gel/snv_hg19_filter_pass"
 
# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"
 
# Loop through each VCF file in the directory
for vcf in "$VCF_DIR"/*.vcf; do
    # Define output file name based on input file
    output_file="$OUTPUT_DIR/$(basename "$vcf" .vcf)_pass.vcf"
 
    # Filter for 'PASS' in the FILTER column (6th column in VCF format)
    # and include header lines starting with ##
    grep -E '^##|PASS' "$vcf" > "$output_file"
done

