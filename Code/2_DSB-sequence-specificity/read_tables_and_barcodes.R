get_barcode_matches <- function(pattern_length, bc_file, max_mismatch=1) {
  # read barcode to sample mappings
  barcode_tab = read.table(bc_file, sep="\t", header=T, stringsAsFactors = F)
  tmp = strsplit(barcode_tab$Sample_name, "_")
  barcode_tab$condition = unlist(sapply(tmp,`[`,2))
  barcode_tab$biosample = unlist(sapply(tmp,`[`,1))
  barcode_tab$condition_short = gsub("Ecoli","", barcode_tab$condition)
  barcode_tab$lib_cond = paste0(barcode_tab$NGS_library,"_", barcode_tab$condition_short)
  
  barcode_tab_wt_mut_ctrl = subset(barcode_tab, biosample=="Caco2" & condition %in% c("EcoliWT","EcoliMut","Ctrl") )
  rownames(barcode_tab_wt_mut_ctrl) = barcode_tab_wt_mut_ctrl$lib_cond
  # only barcodes from conditions used in this analysis
  wt_mut_ctrl_barcodes = unique(barcode_tab_wt_mut_ctrl$Barcode_sequence)
  
  # match patterns to barcodes
  all_patterns = mkAllStrings(c("A","C","G","T"), pattern_length)
  mm = max_mismatch
  pattern_pdict = PDict(DNAStringSet(all_patterns), max.mismatch = mm)
  
  matchfun <- function(x) countPDict(pattern_pdict, DNAString(x), max.mismatch = mm)
  tmp = Map(matchfun, wt_mut_ctrl_barcodes)
  
  bc_forward_matches = lapply(tmp, function(x) all_patterns[which(x>0)])
  #bc_forward_matching_pattern = as.character(all_patterns)[apply(tmp>0, 1, sum) > 0]
  
  matchfun <- function(x) countPDict(pattern_pdict, reverseComplement(DNAString(x)), max.mismatch = mm)
  #tmp = do.call(cbind, Map(matchfun, wt_mut_ctrl_barcodes))
  #bc_revcomp_matching_pattern = as.character(all_patterns)[apply(tmp>0, 1, sum) > 0]
  tmp = Map(matchfun, wt_mut_ctrl_barcodes)
  bc_revcomp_matches = lapply(tmp, function(x) all_patterns[which(x>0)])
  
  return(list("bc_forward_matches"=bc_forward_matches, "bc_revcomp_matches" = bc_revcomp_matches, "barcode_tab_wt_mut_ctrl"=barcode_tab_wt_mut_ctrl))
}

##########################################################################################################
read_enrichment_tables <- function(input_folder, bc_info) {

  # unpack data provided by get_barcode_matches before
  for (n in names(bc_info)) assign(n, bc_info[[n]])
  
  tmp_env = new.env()
  files = list.files(path=input_folder, pattern="*.RData", full.names = T)
  
  all_tables = list()
  for (ff in files) {
    sn = sapply(strsplit(basename(ff), "_"),`[`,2)
    tmp =  unlist(regmatches(sn, regexec("(BB\\d+)(\\S+)", sn)))
    rep_id = tmp[2]
    condition = tmp[3]
    load(ff, envir = tmp_env)
    tmp_tab = as.data.table(get("nt_content", tmp_env))
    tmp_tab$sampleID=sn
    tmp_tab$replicate=rep_id
    tmp_tab$condition = condition
    tmp_1 = strsplit(tmp_tab$condition, "vs")
    tmp_tab$target_condition = unlist(sapply(tmp_1, `[`,1))
    tmp_tab$control_condition = unlist(sapply(tmp_1, `[`,2))
    tmp_tab$bc_target = barcode_tab_wt_mut_ctrl[paste0(tmp_tab$replicate, "_", tmp_tab$target_condition), "Barcode_sequence"]
    #tmp_tab$bc_control = barcode_tab_wt_mut_ctrl[paste0(tmp_tab$replicate, "_", tmp_tab$control_condition), "Barcode_sequence"]
    
    tmp_tab$barcode_match_library = ifelse(as.character(tmp_tab$Nt) %in% bc_forward_matches[[tmp_tab$bc_target[1]]], "BC_forward",ifelse(as.character(tmp_tab$Nt) %in% bc_revcomp_matches[[tmp_tab$bc_target[1]]], "BC_reverse", "no BC"))
    
    tmp_tab$barcode_match_any = ifelse(as.character(tmp_tab$Nt) %in% unlist(bc_forward_matches), "BC_forward", ifelse(as.character(tmp_tab$Nt) %in% unlist(bc_revcomp_matches), "BC_reverse","no BC"))

    # 0/0 -> NaN, but 0 is more sensible
    tmp_tab[, ratio:=ifelse(is.nan(ratio), 0, ratio)]
    # Error propagation for pattern proportions to log2ratios
    # proportions of each motif are approximately poisson distributed. This should be valid for 5nt pattern upwards
    tmp_tab[, sd_r_background:=sqrt(1/size_background * ratio_background * (1-ratio_background))]
    tmp_tab[, sd_r_sample:=sqrt(1/size_sample * ratio_sample * (1-ratio_sample))]
    # SD of ratio
    tmp_tab[, sd_ratio:= sqrt( (sd_r_sample/ratio_sample)**2 + (sd_r_background/ratio_background)**2 ) ]
    # SD of log2ratio
    const = log2(exp(1))
    tmp_tab[, sd_log2ratio:=const*sd_ratio/ratio]
    
    max_logRatio = max(subset(tmp_tab, barcode_match_library=="no BC" & !is.na(logRatio) & !is.infinite(logRatio))$logRatio, na.rm=T)
    tmp_tab[, logRatio_scaled := logRatio / max_logRatio]
    tmp_tab[, logRatio_scale_factor := logRatio / logRatio_scaled] 
    
    # SD of log2ratio_scaled

    tmp_tab[, sd_log2ratio_scaled:= sd_log2ratio / logRatio_scale_factor]

    tmp_tab$rank = nrow(tmp_tab) - rank(tmp_tab$logRatio) + 1
    
    all_tables[[sn]] = tmp_tab
  }
  
  all_nuc_enrichment_tabs = as.data.table(do.call(rbind, all_tables))
  all_nuc_enrichment_tabs[, at_cnt:=sapply(strsplit(as.character(Nt),""), function(x) sum(x %in% c("A","T") ))]
  
  return(all_nuc_enrichment_tabs)
}


