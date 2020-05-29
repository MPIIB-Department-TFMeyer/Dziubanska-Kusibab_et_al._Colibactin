#!/data_genome1/SharedSoftware/R/3.4/lib64/R/bin/exec/R
library('knitr')
library(rmarkdown)

# This script depends on FASTA files containing sequences around DSB extracted before and stored in ../../Data/FASTA_around_DSB

setwd("Compute_enrichment")
source("Run_enrichment_window7_pentanucleotides.R")
source("Run_enrichment_window7_hexanucleotides.R")
source("Run_enrichment_window7_heptanucleotides.R")
source("Run_enrichment_window7_nonanucleotides1.R")
source("Run_enrichment_window7_nonanucleotides2.R")
source("Run_enrichment_window7_nonanucleotides3.R")
source("Run_enrichment_window7_nonanucleotides4.R")
setwd("..")

rmarkdown::render('HN_enrichment_all_replicates.Rmd')
rmarkdown::render('PN_enrichment_all_replicates.Rmd')
rmarkdown::render('AT_content_and_NN_motif_enrichment.Rmd')
rmarkdown::render('HN_distribution_around_nucleosomes.Rmd')
rmarkdown::render('Script_for_scaled_logRatio_HN_table.Rmd')

quit('no')
