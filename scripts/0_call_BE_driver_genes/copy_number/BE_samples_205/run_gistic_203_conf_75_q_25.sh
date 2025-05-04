#!/bin/bash

#! Give your job a name
#SBATCH -J 205_75_gistic
#! How many cores per task?
#SBATCH --cpus-per-task=1
#! How much memory do you need?
#SBATCH --mem=8G

#SBATCH --output=outfile_gistic.out
#SBATCH --error=errorfile_gistic.err

#! How much wallclock time will be required?
#SBATCH --time=10:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=Leanne.Wu@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! General partition
#SBATCH -p general

segfile=/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/data/aks_trio_ross_noora/cna/203_samples_cna_excludeGM_2_exclude_abnormal_cna.seg
#markersfile=`pwd`/examplefiles/markersfile.txt
refgenefile=/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/scripts/gistic/refgenefiles/hg19.mat
#alf=`pwd`/examplefiles/arraylistfile.txt
#cnvfile=`pwd`/examplefiles/cnvfile.txt
output_dir=/scratchb/stlab-icgc/users/wu04/project/bo_gene_list/results/samples_205_excludeGM/conf_75_q_25_no_ploidy_adjust_203_samples

mkdir -p $output_dir

./gistic2 -b $output_dir  \
-seg ${segfile}         \
-refgene ${refgenefile} \
-genegistic 1           \
-smallmem 1             \
-broad 1                \
-armpeel 1              \
-conf 0.75		\
-savegene 1             \
-qvt 0.25		\
1> $output_dir/gistic.log 2>&1
