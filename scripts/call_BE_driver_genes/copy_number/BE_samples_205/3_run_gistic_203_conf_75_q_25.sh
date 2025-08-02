#!/bin/bash
# Script to run GISTIC2 on copy number data for BE samples
# Ensure the GISTIC2 executable is in your PATH or provide the full path to it
segfile=./Genomic_analysis/data/0_call_BE_driver_genes/copy_number/BE_samples_205/cna_for_gistic2_205.seg
#markersfile=`pwd`/examplefiles/markersfile.txt
refgenefile=./Genomic_analysis/scripts/0_call_BE_driver_genes/copy_number/refgenefiles/hg19.mat
#alf=`pwd`/examplefiles/arraylistfile.txt
#cnvfile=`pwd`/examplefiles/cnvfile.txt
output_dir=./Genomic_analysis/results/0_call_BE_driver_genes/copy_number/BE_samples_205/conf_75_q_25_205

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
