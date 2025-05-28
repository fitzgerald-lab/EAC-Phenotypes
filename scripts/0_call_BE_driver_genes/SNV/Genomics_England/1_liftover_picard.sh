#!/bin/bash

#! Give your job a name
#SBATCH -J liftover
#! How many cores per task?
#SBATCH --cpus-per-task=1
#! How much memory do you need?
#SBATCH --mem=8G

#SBATCH --output=outfile.out
#SBATCH --error=errorfile.err

#! How much wallclock time will be required?
#SBATCH --time=05:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=Leanne.Wu@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! General partition
#SBATCH -p general

VCF_DIR="/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/gel/snv"

# Define output directory for lifted VCFs
OUTPUT_DIR="/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/gel/snv_hg19"

# Location of the Picard jar
PICARD_JAR="/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/scripts/picard/picard.jar"

# Chain file for liftover
CHAIN_FILE="/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/scripts/picard/hg38ToHg19.over.chain"

# Reference sequence fasta
REFERENCE_SEQUENCE="/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/scripts/picard/hg19.fa"


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
