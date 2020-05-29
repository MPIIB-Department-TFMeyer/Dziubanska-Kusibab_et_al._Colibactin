#!/data_genome1/SharedSoftware/R/3.4/lib64/R/bin/exec/R
library(knitr)
library(rmarkdown)

# This script depends on FASTA files containing sequences around DSB extracted before and stored in ../../Data/FASTA_around_DSB

setwd("AsiSI_validation")
rmarkdown::render('import_break_ends.Rmd')
rmarkdown::render('Break_end_distribution_in_AsiSI_motif.Rmd')
setwd("..")

rmarkdown::render('Break_ends_analysis_in_selected_motifs.Rmd')
rmarkdown::render('Frequencies_of_nucleotide_changes_CRC_WGS_in_CDM.Rmd')

quit('no')
