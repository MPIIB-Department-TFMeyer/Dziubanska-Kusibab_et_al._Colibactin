# Overview
Scripts and data for analyses of Colibactin associated DNA damage in Dziubanska-Kusibab &amp; Berger et al. (DOI: 10.1038/s41591-020-0908-2 )

# Requirements 

## Hardware and Software

All analyses were run on a 48 core server with 256 Gbytes of RAM, running Ubuntu 16.04. Data, code and results require at least approximately 50G of available hard disk space. 

To reproduce the analysis you will need some or all of the following software: 

  - R (>= v3.4.0) [cran.r-project.org]
  - MEME suite (>= v5.0.5)    [http://meme-suite.org/doc/download.html]

The following R packages should be installed from either CRAN or Bioconductor: 

  - BiocParallel
  - Biostrings
  - BSgenome.Hsapiens.UCSC.hg19
  - BSgenome.Hsapiens.UCSC.hg38
  - data.table
  - deconstructSigs
  - directlabels
  - DNAshapeR
  - dplyr
  - GenomicFeatures
  - GenomicRanges
  - ggplot2
  - ggpubr
  - ggseqlogo
  - ggthemes
  - grid
  - gridExtra
  - Hmisc
  - knitr
  - naturalsort
  - pheatmap
  - randomForest
  - readxl
  - reshape2
  - rgl
  - rmarkdown
  - robust
  - robustbase
  - scales
  - seqRFLP
  - stringr
  - VariantAnnotation

## External data

Some data used in this analysis has not been included in this repository due to license limitations or size. All external data sets are listed under ./Data/External, accompanied by   source URLs and further scripts and instructions for processing. 

In addition, for generation of DSB-contect FASTA files a human genome sequence is required (we used human_g1k_v37_decoy.fasta from the GATK bundle https://gatk.broadinstitute.org/hc/en-us, https://github.com/broadinstitute/gatk/releases).

# Run the analysis

Each subfolder under ./Code contains a file *run_all.sh* that should run all scripts of a given part of the study and will produce result tables and figures under /Results. 

# Figures

Mappings of figures in the published manuscript to individual scripts can be found in file *figure_mapping.txt*

