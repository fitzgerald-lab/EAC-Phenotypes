This repository contains the data files and analysis code used in **"Integrated epidemiological and molecular data yields insights into the relationship between precancerous Barrett’s and oesophageal adenocarcinoma"** published in [journal name] on [date]. 

# Overview
The files are organized into three main folders:

- data: contains de-identified example input data.
- scripts: contains analysis code to reproduce the analysis.
- results: filtered results generated from the soruce files.

# Set up
## Tools and libraries

| Tool & Version              | Purpose                      | Input data          | Output results                   | Source                                    |
|:----------------------------|:-----------------------------|:--------------------|:---------------------------------|:------------------------------------------|
| Strelka (2.9.4 for Genomics England WGS data, 2.0.15 for other WGS data) | Call somatic mutations and indels  | Aligned sequencing data (.bam files) | Somatic SNVs and indels          | https://github.com/Illumina/strelka       |
| ASCAT 2.1                   | Call copy number for data not from Genomics England | Aligned .bam files  | Allele-specific copy number profiles     | https://github.com/VanLoo-lab/ascat       |
| Manta                       | Call structural variants     | Aligned .bam files  | Structural variant calls (SVs)   | https://github.com/Illumina/manta         |
| Canvas 1.38.0.1554          | Infer copy number for Genomics England data | Aligned .bam files  | Allele-specific copy number profiles | https://github.com/Illumina/canvas        |
| Picard 3.1.1 (LiftoverVcf)  | Convert VCF files from hg38 to hg19 coordinates      | VCF files (hg38)    | Lifted VCF files (hg19)          | https://broadinstitute.github.io/picard/  |
| GISTIC2.0                   | Identify recurrent amplifications/deletions across the Barrett’s cohort | Copy number segmentation files | G-scores, peak regions, significant genes | https://github.com/broadinstitute/gistic2 |
| dNdScv 0.1.0                | Identify positively selected | SNV/indel mutation  | List of significant driver genes | https://github.com/im3sanger/dndscv       |
| SigProfilerExtractor & deconstructSigs | Identify and quantify mutational signatures | VCF files with SNVs and  indels | Sample-specific mutational signatures and contributions | https://github.com/AlexandrovLab/SigProfilerExtractor; https://github.com/raerose01/deconstructSigs |
| Amplicon Architect 1.2 &  Amplicon Classifier 0.4.13 | Reconstruct and classify amplified genomic regions (e.g., ecDNA, BFB) | Aligned .bam files  | Amplicon structures and classifications | https://github.com/virajbdeshpande/AmpliconArchitect; https://github.com/AmpliconSuite/AmpliconClassifier |


Required libraries (any version)
```
library(tidyverse)
library(finalfit)
library(broom)
library(performance)
library(table1)
library(kableExtra)
```
Local computing environment

R version 4.3.2, Platform: x86_64-conda-linux-gnu (64-bit), macOS Ventura 13.3.1
