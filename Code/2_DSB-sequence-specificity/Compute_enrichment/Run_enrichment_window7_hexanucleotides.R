rm(list=ls())

########################################## read fasta & convert to df ######################

source("./compute_enrichment.R")

############################################################################################

input_folder = "../../../Data/FASTA_around_DSB/dsb_pm_7/"
output_folder = "../../../Data/Sequence_enrichment/Hexanucleotides/"
if ( ! file.exists(output_folder) ) dir.create(output_folder)

oligo_size = 6

# BB76&BB77

compute_enrichment(testfile=file.path(input_folder,"BB77_WT_window7.fasta"),
                   controlfile = file.path(input_folder,"BB77_Mut_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB77WTvsMut", 
                   output_folder = output_folder)

compute_enrichment(testfile=file.path(input_folder,"BB77_WT_window7.fasta"), 
                   controlfile =file.path(input_folder,"BB76_Ctrl_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB77WTvsCtrl", 
                   output_folder =output_folder)

compute_enrichment(testfile=file.path(input_folder,"BB77_Mut_window7.fasta"), 
                   controlfile =file.path(input_folder,"BB76_Ctrl_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB77MutvsCtrl", 
                   output_folder =)

# BB78vsBB79

compute_enrichment(testfile=file.path(input_folder,"BB79_WT_window7.fasta"), 
                   controlfile = file.path(input_folder,"BB79_Mut_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB79WTvsMut", 
                   output_folder = output_folder)

compute_enrichment(testfile=file.path(input_folder,"BB79_WT_window7.fasta"), 
                   controlfile =file.path(input_folder,"BB78_Ctrl_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB79WTvsCtrl", 
                   output_folder =output_folder)

compute_enrichment(testfile=file.path(input_folder,"BB79_Mut_window7.fasta"), 
                   controlfile =file.path(input_folder,"BB78_Ctrl_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB79MutvsCtrl", 
                   output_folder =output_folder)


# BB82&BB83

compute_enrichment(testfile=file.path(input_folder,"BB83_WT_window7.fasta"), 
                   controlfile = file.path(input_folder,"BB83_Mut_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB83WTvsMut", 
                   output_folder = output_folder)

compute_enrichment(testfile=file.path(input_folder,"BB83_WT_window7.fasta"), 
                   controlfile =file.path(input_folder,"BB82_Ctrl_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB83WTvsCtrl", 
                   output_folder =output_folder)

compute_enrichment(testfile=file.path(input_folder,"BB83_Mut_window7.fasta"), 
                   controlfile =file.path(input_folder,"BB82_Ctrl_window7.fasta"), 
                   oligo_size=oligo_size, 
                   outputname = "Enrichment_BB83MutvsCtrl", 
                   output_folder =output_folder)

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
