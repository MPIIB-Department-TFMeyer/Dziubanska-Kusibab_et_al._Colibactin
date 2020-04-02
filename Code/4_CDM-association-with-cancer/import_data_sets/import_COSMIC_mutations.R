#setwd("/path/to/this/folder")
library(VariantAnnotation)
library(data.table)

# This gives the genomic coordinates and REF/ALT sequences for unique mutations
cosmic_sites = readVcf("../../../Data/External/COSMIC_mutations/CosmicCodingMuts.vcf.gz")

cosmic_ranges = rowRanges(cosmic_sites)
cosmic_ranges$CNT = info(cosmic_sites)[["CNT"]]
cosmic_ranges$GENE = info(cosmic_sites)[["GENE"]]
cosmic_ranges$STRAND = info(cosmic_sites)[["STRAND"]]
cosmic_ranges$CDS = info(cosmic_sites)[["CDS"]]
cosmic_ranges$AA = info(cosmic_sites)[["AA"]]
cosmic_ranges$SNP = info(cosmic_sites)[["SNP"]]

rl = nchar(cosmic_ranges$REF)
al = nchar(unstrsplit(CharacterList(cosmic_ranges$ALT), ","))

cosmic_ranges$variant_type = ifelse(rl!=1 | al!=1, "InDel", "SNV")

# most of the targeted mutation data has no genomic position in the DB
# only a fraction (200k) of them match the variants in the vcf
#cosmic_targeted = fread("zcat CosmicCompleteTargetedScreensMutantExport.tsv.gz")

# this table matches > 95% of variant sites in the vcf
# it contains individual mutations *and* for individual samples, but not the genomic changes and annotations in an optimal format
cosmic_ngs = fread("zcat ../../../Data/External/COSMIC_mutations/CosmicGenomeScreensMutantExport.tsv.gz")
setkey(cosmic_ngs, "Mutation ID")

shared_ids = unique(cosmic_ngs$`Mutation ID`[cosmic_ngs$`Mutation ID` %in% names(cosmic_ranges)])

tmp1 = subset(cosmic_ngs, cosmic_ngs$`Mutation ID` %in% shared_ids)
cosmic_mutations = cosmic_ranges[tmp1$`Mutation ID`]
mcols(cosmic_mutations) = cbind(mcols(cosmic_mutations), tmp1)
cosmic_mutations$stype = with(cosmic_mutations, paste(Primary.site, Primary.histology,Histology.subtype.1, sep="_"))

tmp = fread("../../../Data/COSMIC_mutations/COSMIC_subtype_to_acronym.csv", sep="\t")
setkey(tmp, "type")
tmp[tmp$acronym==""]$acronym = NA

cosmic_mutations$entity_drivers=rep(NA_character_, length(cosmic_mutations))
cosmic_mutations$entity_drivers = tmp[as.character(cosmic_mutations$stype)]$acronym

tmp_pole_flag = which(cosmic_mutations$GENE=="POLE" & !cosmic_mutations$Mutation.Description=="Substitution - coding silent")
pole_mutated_cases_cosmic = unique(cosmic_mutations[tmp_pole_flag]$ID_sample)
tmp_pold_flag = which(cosmic_mutations$GENE=="POLD1" & !cosmic_mutations$Mutation.Description=="Substitution - coding silent")
pold_mutated_cases_cosmic = unique(cosmic_mutations[tmp_pold_flag]$ID_sample)


cosmic_mutations$POLE_mutated = ifelse(cosmic_mutations$ID_sample %in% pole_mutated_cases_cosmic, T, F)
cosmic_mutations$POLD_mutated = ifelse(cosmic_mutations$ID_sample %in% pold_mutated_cases_cosmic, T, F)

cosmic_mutations$TCGA = grepl("TCGA", cosmic_mutations$`Sample.name`)

save(cosmic_mutations, file="../../../Data/COSMIC_mutations/COSMIC_mutation_ranges.Rdata")
