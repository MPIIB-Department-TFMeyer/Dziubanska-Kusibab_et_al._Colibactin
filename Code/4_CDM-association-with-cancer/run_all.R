#!/data_genome1/SharedSoftware/R/3.4/lib64/R/bin/exec/R
library(knitr)
library(rmarkdown)

# This script depends on FASTA files containing sequences around DSB extracted before and stored in ../../Data/FASTA_around_DSB
setwd("import_data_sets")
source("import_som_mutations_Giannakis_WXS.R")
source("import_MC3_TCGA_mutations.R")
source("import_COSMIC_mutations.R")
setwd("..")

setwd("Pattern_positions_genomic")
rmarkdown::render('Generate_HN_Positions_4096.Rmd')
rmarkdown::render('Generate_HN_Positions_hg38.Rmd')
rmarkdown::render('Generate_PN_Positions_hg38.Rmd')
setwd("..")

setwd("Giannakis_CDM_enrichment")
rmarkdown::render('Pattern_enrichment_around_variants_Giannakis_COREAD_quartiles_and_HM_HN_4096_v2.Rmd')
rmarkdown::render('Pattern_enrichment_around_variants_Giannakis_COREAD_quartile_and_HM_HN_4096_v2_Plots_final.Rmd')
setwd("..")

setwd("SBS_Derivation")
rmarkdown::render('Theoretical_AAWWTT_SBS_signature.Rmd')
setwd("..")

setwd("MC3_CDM_enrichment/")
rmarkdown::render('Pattern_enrichment_around_variants_MC3_PanCancer_quartiles_and_HM_HN_4096_v2.Rmd')
rmarkdown::render('Pattern_enrichment_around_variants_MC3_PanCancer_quartile_and_HM_HN_4096_v2_Plots.Rmd')
setwd("..")

setwd("CRC_WGS_CDM_enrichment")
rmarkdown::render('Finnish_WGS_cohort_HN.Rmd')
rmarkdown::render('Finnish_WGS_cohort_PN.Rmd')
rmarkdown::render('Pattern_enrichment_around_variants_Finnish_COREAD_quartile_and_HM_HN_64_Plots.Rmd')
rmarkdown::render('Pattern_enrichment_around_variants_Finnish_COREAD_quartile_and_HM_PN_Plots.Rmd')
rmarkdown::render('WGS_cohort_SBS_signature_analysis_COSMICv3_and_LeeSix.Rmd')
rmarkdown::render('Finnish_cohort_further_correlations.Rmd')
setwd("..")

setwd("Driver_genes/")
rmarkdown::render('Driver_genes.Rmd')
rmarkdown::render('Driver_gene_hits_PanCancer_Drivers_PN_and_HN_COSMIC.Rmd')
rmarkdown::render('Driver_gene_hits_PanCancer_Drivers_PN_and_HN_MC3.Rmd')
rmarkdown::render('Analysis_of_hits_at_HN_and_PN_in_COSMIC_v2.Rmd')
rmarkdown::render('Analysis_of_hits_at_HN_and_PN_in_MC3_v2.Rmd')
setwd("..")

setwd("Signature_contamination_analysis")
rmarkdown::render('PCAWG_WGS_cohort_signature_analysis_global_Alexandrov_2018.Rmd')
rmarkdown::render('PentanucleotideChanges.Rmd')
setwd("..")

quit('no')
