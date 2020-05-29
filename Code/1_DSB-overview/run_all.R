#!/data_genome1/SharedSoftware/R/3.4/lib64/R/bin/exec/R
library('knitr')
library(rmarkdown)

# Please make sure that the MEME suite (>= v5.05) is installed and that paths in the DREME scripts have been adjusted

system("pushd DREME_E.coli; bash SCRIPT_RUN_DREME; popd")
setwd("DREME_E.coli/analysis")
system("bash run_R_DREME_pos_analysis.sh")

rmarkdown::render('DREME_motif_overview.Rmd')
rmarkdown::render('DREME_positional_analysis.Rmd')
setwd("../..")

system("pushd DREME_Etoposide; bash SCRIPT_RUN_DREME; popd")
setwd("DREME_Etoposide/analysis")
system("bash run_R_DREME_pos_analysis.sh")

rmarkdown::render('DREME_motif_overview.Rmd')
rmarkdown::render('DREME_positional_analysis.Rmd')
setwd("../..")

quit('no')


