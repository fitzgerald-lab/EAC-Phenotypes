This repository contains the data files and analysis code used in **"Integrated epidemiological and molecular data yields insights into the relationship between precancerous Barrett’s and oesophageal adenocarcinoma"** published in [journal name] on [date]. 

# Overview
The files are organized into three main folders:

- data: contains de-identified example input data.
- scripts: contains analysis code to reproduce the analysis.
- results: filtered results generated from the soruce files.

# Data availability

| Data              | Description                  | Source(s)              | Demo data              |
|:------------------|:-----------------------------|:-----------------------|:-----------------------|
| Whole-genome and -exome sequencing data | Genomic data inlcuded in the current study. Access request may be required per repository | https://ega-archive.org/datasets/EGAD00001015467; https://ega-archive.org/studies/EGAS00001003702; https://research.genomicsengland.co.uk/research-registry/browse/ | [↗️]([https://github.com/fitzgerald-lab/EAC-Phenotypes/tree/main/data/genomic]) |
| Clinical data | Clinical and epidemiological data for regression models | Full dataset request via https://www.occams.org.uk/index.html | [↗️](https://github.com/fitzgerald-lab/EAC-Phenotypes/tree/main/data/epidemiological)



# Set up
## Tools and libraries

| Tool & Version              | Purpose                      | Input data          | Output results                   | Source(s)                                    |
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
| ShatterSeek 1.1      | To identify chromothripsis events    | Structural variants + copy number profiles  | Predicted chromothripsis events | https://github.com/parklab/ShatterSeek   |
| GATK Mutect2 4.1.7.0 | Call somatic mutations in WES samples    | Aligned .bam files      | Somatic SNV and indel calls   | https://gatk.broadinstitute.org/hc/en-us/articles/360036485152-Mutect2 |
| ClonEvol 0.99.11     | Construct phylogenetic trees of tumour evolution  | Subclonal populations data   | Reconstructed clonal phylogenies | https://github.com/hdng/clonevol         |
| PyClone 0.13.1       | Infer subclonal populations from SNV data   | Somatic mutation table  | Subclonal clusters with VAF and CCF estimates  | https://github.com/Roth-Lab/pyclone      |
| REVOLVER 1.0.0       | Infer evolutionary trajectories across patient cohort | Subclonal mutation profiles from multiple regions  | Evolutionary trees, clusters, and recurrent trajectories | https://github.com/caravagnalab/revolve  |
| Space Ranger 3.1.3   | Process 10x Genomics Visium HD spatial transcriptomics data | Raw FASTQ and image data | Aligned gene expression matrices and spatial metadata | http://github.com/10XGenomics/spaceranger    |
| Loupe Browser 8.1.2  | Visualize and analyze spatial transcriptomics data | Space Ranger output files     | Interactive spatial gene expression plots  | https://www.10xgenomics.com/support/software/loupe-browser/latest     |
| Scanpy 1.11.1        | Analyze spatial gene expression matrices          | Gene expression matrix        | Cluster annotations and spatial plots | https://github.com/scverse/scanpy |

Required R libraries (any version)

```{r}
# data wrangling
library(tidyverse) 

# imputation
library(mice)

# model building
library(finalfit)
library(broom)
library(performance)

# organizing results
library(table1)
library(kableExtra)
```

## Local computing environment

R version 4.3.2, Platform: x86_64-conda-linux-gnu (64-bit), macOS Ventura 13.3.1
