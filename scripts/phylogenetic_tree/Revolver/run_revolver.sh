#!/bin/bash

#! Give your job a name
#SBATCH -J revolver 
#! How many cores per task?
#SBATCH --cpus-per-task=2
#! How much memory do you need?
#SBATCH --mem=64G

#SBATCH --output=outfile_jack.out
#SBATCH --error=errorfile_jack.err

#! How much wallclock time will be required?
#SBATCH --time=20:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=Leanne.Wu@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! General partition
#SBATCH -p epyc


Rscript /mnt/scratchc/stlab-icgc/users/wu04/project/bo_gene_list/scripts/phylogenetic_tree/run_revolver.R
