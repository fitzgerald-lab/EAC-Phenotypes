#!/bin/bash

#! Give your job a name
#SBATCH -J indel_dNdS 
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

Rscript /scratchb/stlab-icgc/users/wu04/project/bo_gene_list/scripts/samples_211/from_vep_indel_to_dNdS.R
