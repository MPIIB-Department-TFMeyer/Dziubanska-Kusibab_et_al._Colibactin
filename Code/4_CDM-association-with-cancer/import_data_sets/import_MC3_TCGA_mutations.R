#setwd("/path/to/this/folder")

library(data.table)
library(GenomicRanges)
library(GenomicFeatures)

data = fread("../../../Data/External/MC3_TCGA/MC3_variants_with_trinuc_context.txt")

# Check for donors with more than one tumor sample
ed_tmp = data[, .(bc_count = length(unique(Tumor_Sample_Barcode)), cancer_entity=unique(cancer_entity)), by="TCGA_donor"]
#  as.data.frame(unique(data[, c("TCGA_donor", "cancer_entity")]))
setkey(ed_tmp, "TCGA_donor")
# There are about 70 cases across several cancer entitiies
sample_with_gt1_bc = ed_tmp[ed_tmp$bc_count>1]$TCGA_donor

# Check for donors with more than one normal sample - in those cases variant calls are effectively duplicated
ed_tmp = data[, .(bc_count = length(unique(Matched_Norm_Sample_Barcode)), cancer_entity=unique(cancer_entity)), by="TCGA_donor"]
#  as.data.frame(unique(data[, c("TCGA_donor", "cancer_entity")]))
sample_with_gt1_control = ed_tmp[ed_tmp$bc_count>1]$TCGA_donor
setkey(ed_tmp, "TCGA_donor")

# there are about 490 cases with > 1 normal control
# In all cases there is one from blood, one from normal tissue
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
# 10 	Blood Derived Normal 	NB
# 11 	Solid Tissue Normal 	NT

tmp = data[, .(count=.N), by=c("TCGA_donor","Matched_Norm_Sample_Barcode")]
tmp2 = strsplit(tmp$Matched_Norm_Sample_Barcode,"-")
tmp3 = unlist(sapply(tmp2, `[`,4))
tmp[, control_type:= substr(tmp3,1,2) ]
tmp[, count:=.N, by="TCGA_donor"]

tmp[, selected_control:=ifelse(count==1, TRUE, ifelse(control_type==10, TRUE, FALSE))]
setkey(tmp, "TCGA_donor")
selected_germline_controls = tmp[tmp$selected_control]

# Now remove all duplicate calls due to >1 normal controls
data_orig = data
data = data[data$Matched_Norm_Sample_Barcode == selected_germline_controls[data$TCGA_donor]$Matched_Norm_Sample_Barcode]

# Trinuc changes annotation
trinuc_change_anno = as.data.frame(unique(data[data$Variant_Type=="SNP", c("trinuc_context", "trinuc_change", "nuc_change")]))
rownames(trinuc_change_anno) = trinuc_change_anno$trinuc_change

# all cacner entities/subtypes
mc3_mutation_entities = unique(sort(data$cancer_entity)) 

# Sample exclusion list from Bailey et al. 2018
sample_exclusion_list = fread("../../../Data/External/MC3_TCGA/TCGA_exlcuded_cases_Bailey_et_al_2018.txt", sep="\t")

sample_exclusion_list_final = subset(sample_exclusion_list, sample_exclusion_list$`Reason removed`!="RNA degradation")
setkey(sample_exclusion_list_final, "Sample barcode")
#data[, sample_exclusion:=rep(NA_character_, nrow(data))]
data[, sample_exclusion:=ifelse(Tumor_Sample_Barcode %in% sample_exclusion_list_final$`Sample barcode`, sample_exclusion_list_final[Tumor_Sample_Barcode]$`Reason removed`, NA)]

# POLE/POLD1 mutated cases show often an hypermutator phenotype
pole_mutated_cases = subset(data, Hugo_Symbol=="POLE" & IMPACT %in% c("MODERATE", "HIGH"))$TCGA_donor
pold_mutated_cases = subset(data, Hugo_Symbol=="POLD1" & IMPACT %in% c("MODERATE", "HIGH"))$TCGA_donor

# Keep only sample with a single tumor sample and SNV/InDel mutations
mc3_data_snp_indel = data[data$Variant_Type %in% c("SNP", "INS","DEL") & !TCGA_donor %in% sample_with_gt1_bc & !cancer_entity==""]

# Genomic ranges objects
mc3_mut_ranges_snp_indel = with(mc3_data_snp_indel, GRanges(Chromosome, IRanges(Start_Position, End_Position), strand="*", Reference_Allele, Tumor_Seq_Allele2, nuc_change, donorID = factor(TCGA_donor), trinuc_change,  mutation_type = Variant_Type, cancer_entity = cancer_entity ) )
mc3_mut_ranges_snp_indel$var_id = with(mc3_data_snp_indel, paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, cancer_entity, sep="_"))

# Exclude non-exonic regions 
gencode28 = loadDb("../../../Data/External/Gencode/GENCODEv28lift37_basic_hg19_TxDB.db")

ee = exonsBy(gencode28, "gene")
sc = standardChromosomes(ee)
sc = sc[sc!="MT"]
ee = keepSeqlevels(ee, sc, pruning.mode = "coarse")

exon_ranges = reduce(ee, ignore.strand=TRUE)

oo = findOverlaps(mc3_mut_ranges_snp_indel, exon_ranges, ignore.strand=T)

mc3_mut_ranges_snp_indel_exonic_only = mc3_mut_ranges_snp_indel[unique(queryHits(oo))]

# sample annotation including sample exclusion by Bailey et al and POLE/POLD1 mutation status
ed_MC3 = as.data.frame(unique(mc3_data_snp_indel[, c("TCGA_donor", "cancer_entity", "sample_exclusion")]))
rownames(ed_MC3) = ed_MC3$TCGA_donor
ed_MC3$POLE_mutation = ifelse(ed_MC3$TCGA_donor %in% pole_mutated_cases, T, F)
ed_MC3$POLD_mutation = ifelse(ed_MC3$TCGA_donor %in% pold_mutated_cases, T, F)

#table(ifelse(is.na(ed_final$sample_exclusion),"n.a.", ed_final$sample_exclusion), ed_final$POLD_mutation|ed_final$POLE_mutation, dnn=c("Exclusion Bailey et al.", "POLE/POLD1 mutation"))

total_mut_counts_ts = mc3_data_snp_indel[, .(count = .N), by = c("TCGA_donor", "Variant_Type")]
total_mut_counts_tab = dcast.data.table(total_mut_counts_ts, TCGA_donor ~ Variant_Type, value.var="count")
total_mut_counts_tab = total_mut_counts_tab[,lapply(.SD,function(x){ifelse(is.na(x),0,x)})]
total_mut_counts_tab[, INDEL:=INS+DEL]
total_mut_counts_tab[, total_variants:=SNP+INS+DEL]
setkey(total_mut_counts_tab, "TCGA_donor")

# Add variant counts per sample
ed_MC3 = merge(ed_MC3, total_mut_counts_tab, by="TCGA_donor", all.x=T)

# Count variants still included after restricting to exonic regions
tmp = data.table(TCGA_donor = mc3_mut_ranges_snp_indel_exonic_only$donorID)[, .(count_exonic=.N), by="TCGA_donor"]
setkey(tmp, "TCGA_donor")
ed_MC3 = merge(ed_MC3, tmp, by="TCGA_donor", all.x=T)

fwrite(mc3_data_snp_indel, file="../../../Data/MC3_TCGA_WXS/MC3_SNP_InDel_mutation_table.txt", sep="\t")
save(mc3_mut_ranges_snp_indel_exonic_only, file="../../../Data/MC3_TCGA_WXS/MC3_SNP_InDel_mutation_ranges.Rdata")
save(ed_MC3, trinuc_change_anno, mc3_mutation_entities, file="../../../Data/MC3_TCGA_WXS/MC3_sample_anno.Rdata")
