rm(list=ls())

########################################## read fasta & convert to df ######################

source("./compute_enrichment.R")

############################################################################################

input_folder = "../../../Data/FASTA_around_DSB/dsb_pm_7/"
output_folder = "../../../Data/Sequence_enrichment/Nonanucleotides/"
if ( ! file.exists(output_folder) ) dir.create(output_folder)

oligo_size = 9

# BB191&BB192

compute_enrichment(testfile=file.path(input_folder,"BB192_WT_window7.fasta"), 
                   controlfile = file.path(input_folder,"BB192_Mut_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB192WTvsMut", 
                   output_folder = output_folder)

compute_enrichment(testfile=file.path(input_folder,"BB192_WT_window7.fasta"), 
                   controlfile =file.path(input_folder,"BB191_Ctrl_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB192WTvsCtrl", 
                   output_folder =output_folder)

compute_enrichment(testfile=file.path(input_folder,"BB192_Mut_window7.fasta"), 
                   controlfile =file.path(input_folder,"BB191_Ctrl_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB192MutvsCtrl", 
                   output_folder =output_folder)
