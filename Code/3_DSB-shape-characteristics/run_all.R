#!/data_genome1/SharedSoftware/R/3.4/lib64/R/bin/exec/R
library('knitr')
library(rmarkdown)

# This script depends on FASTA files containing sequences around DSB extracted before and stored in ../../Data/FASTA_around_DSB

rmarkdown::render('Prepare_DNA_Shape_and_Enrichment_tables.Rmd')
rmarkdown::render('Sequence_and_Shape_correlation.Rmd')
setwd("Multivariate")
rmarkdown::render('RF_and multivariate_analysis_6nt_with_stiffness_4reps.Rmd')
# Run 3D_final.R on a computer with X in order to allow rgl to work
setwd("..")
setwd("Complete_model")
rmarkdown::render('Pattern_analysis_shape_influence.Rmd')
setwd("..")
quit('no')
