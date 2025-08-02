#!/bin/bash

# Input =  ../snv_hg19/*.vcf
# Output =  ../snv_hg19_filter_pass/*_pass.vcf

# Directory containing your VCF files
VCF_DIR="/path/to/input/vcf_files"
 
# Output directory for filtered VCF files
OUTPUT_DIR="/path/to/output/filtered_vcf_files"
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

# End of script
