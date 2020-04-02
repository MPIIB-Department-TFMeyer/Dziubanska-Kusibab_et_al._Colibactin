
lee_six_sigs = read.table("/data_genome1/public_data/Stratton_crypt_signatures/git/colon_microbiopsies/signature_extraction/subsitutions_hdp_signature_extraction/Signature_category_counts.txt", sep="\t", header=T)

lee_six_mat = as.matrix(lee_six_sigs[, 2:ncol(lee_six_sigs)])
tmp_c = as.character(lee_six_sigs$X)
nuc_change = gsub("\\.",">", substr(tmp_c, 1, 3))
trinuc_change = paste0(substr(tmp_c,8,8),"[",nuc_change,"]",substr(tmp_c,10,10))

rownames(lee_six_mat) = trinuc_change

lee_six_mat_new_sigs = lee_six_mat[, c(paste0("N",1:5))]

new_sigs = read.csv(file="/data_genome2/public_data/Alexandrov_ICGC/v2018/sigProfiler_SBS_signatures_2018_03_28.csv")
new_sigs_fixed = as.data.frame(t(new_sigs[, 3:ncol(new_sigs)]))
colnames(new_sigs_fixed) = apply(new_sigs[, 1:2], 1, function(x) {paste0(substr(x[2],1,1),"[", x[1],"]", substr(x[2],3,3))})

alexandrov_and_lee_six_sigs = rbind(t(lee_six_mat_new_sigs), new_sigs_fixed[,rownames(lee_six_mat_new_sigs)])


##################################################################
# All trinuc changes
##################################################################
tmp = data.frame(tc = trinuc_change, stringsAsFactors = F)
tmp$nc = substr(tmp$tc, 3,5)
tmp$context = paste0(substr(tmp$tc,1,1), substr(tmp$tc,7,7))
tmp = tmp[order(tmp$nc, tmp$context),]
trinuc_changes_ordered = tmp
rownames(trinuc_changes_ordered) = trinuc_changes_ordered$tc
trinuc_changes_ordered$tn= with(trinuc_changes_ordered, paste0(substr(context,1,1), substr(nc,1,1),substr(context,2,2)))

##################################################################
# expected trinucleotides affected by mutations in colibactin associated motifs
##################################################################
expected_mut_CAM = subset(trinuc_changes_ordered, !grepl("NA",nc))
expected_mut_CAM$CAM5 = 0
# Manually computed frequencies for TN in 5nt + flanking N for AAAAT/AAATT. 10 total trinucleotide positions, 3 possible NC
expected_mut_CAM[expected_mut_CAM$tn=="ATA",]$CAM5 <- 0.25/10/3
expected_mut_CAM[expected_mut_CAM$tn=="ATC",]$CAM5 <- 0.25/10/3
expected_mut_CAM[expected_mut_CAM$tn=="ATG",]$CAM5 <- 0.25/10/3
expected_mut_CAM[expected_mut_CAM$tn=="ATT",]$CAM5 <- 3.25/10/3
expected_mut_CAM[expected_mut_CAM$tn=="TTA",]$CAM5 <- 0.75/10/3
expected_mut_CAM[expected_mut_CAM$tn=="TTC",]$CAM5 <- 0.75/10/3
expected_mut_CAM[expected_mut_CAM$tn=="TTG",]$CAM5 <- 0.75/10/3
expected_mut_CAM[expected_mut_CAM$tn=="TTT",]$CAM5 <- 3.75/10/3
expected_mut_CAM$CAM5_factor = ifelse(expected_mut_CAM$CAM5==0, 1, expected_mut_CAM$CAM5)
expected_mut_CAM$tc = factor(expected_mut_CAM$tc, levels=trinuc_changes_ordered$tc)

expected_mut_CAM$CAM6 = 0
# Manually computed frequencies for TN in 5nt + flanking N for AAAAT/AAATT. 10 total trinucleotide positions, 3 possible NC
expected_mut_CAM[expected_mut_CAM$tn=="ATA",]$CAM6 <- 2/18/3
expected_mut_CAM[expected_mut_CAM$tn=="ATT",]$CAM6 <- 5/18/3
expected_mut_CAM[expected_mut_CAM$tn=="TTA",]$CAM6 <- 1.5/18/3
expected_mut_CAM[expected_mut_CAM$tn=="TTC",]$CAM6 <- 1.5/18/3
expected_mut_CAM[expected_mut_CAM$tn=="TTG",]$CAM6 <- 1.5/18/3
expected_mut_CAM[expected_mut_CAM$tn=="TTT",]$CAM6 <- 6.5/18/3
expected_mut_CAM$CAM6_factor = ifelse(expected_mut_CAM$CAM6==0, 1, expected_mut_CAM$CAM6)

expected_mut_CAM$CAM6_hotspot_factor = 1
# Manually computed frequencies for TN in 5nt + flanking N for AAAAT/AAATT. 10 total trinucleotide positions, 3 possible NC
expected_mut_CAM[expected_mut_CAM$tn=="ATT" & expected_mut_CAM$nc=="T>C",]$CAM6_hotspot_factor <- 4 * 3/4
expected_mut_CAM[expected_mut_CAM$tn=="TTT" & expected_mut_CAM$nc=="T>C",]$CAM6_hotspot_factor <- 4 * 1/4

expected_mut_CAM$CAM6_hotspot = expected_mut_CAM$CAM6 * expected_mut_CAM$CAM6_hotspot_factor
expected_mut_CAM$CAM6_hotspot = expected_mut_CAM$CAM6_hotspot/sum(expected_mut_CAM$CAM6_hotspot)

CAM6_trinucleotides = unique(subset(expected_mut_CAM, CAM6>0)$tn)
