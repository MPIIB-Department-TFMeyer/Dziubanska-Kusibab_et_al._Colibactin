require(Biostrings)
require(GenomicRanges)

# We assume that query and target regions have a direction (i.e. stranded-ness)
# target regions are typically short oligonucleotide motifs
# query regions are estimated double strand break sites derived either from a + or - stranded read and are
# given as single nucleotide positions 5' of the break site
# We assume that the point of reference in target regions is the center of the provided regions
# We also assume that palindromic target regions are marked as such (those will only mapped to the + strand)

compute_distances <- function(query_region, target_regions, maxgap=1, break_between_nucleotides=F, palindrome_ref_strand = "+", use_center_ref_position = F) {
  
  pp = target_regions
  # for motifs with odd number of nucleotides, center is the middle nucleotide. 
  # for motifs with even number of nucleotides, center is the nucleotide following the (n/2)th nucleotide
  # positions are always counted in 5' -> 3' direction
  
  if (use_center_ref_position) {
    pp$region_ref_position = ifelse(strand(pp)=="+", 
                              ifelse((width(pp) %% 2)==0, start(pp) + width(pp)/2, start(pp) + (width(pp)-1)/2), 
                              ifelse((width(pp) %% 2)==0, end(pp) - width(pp)/2,   end(pp) - (width(pp)-1)/2) )
    
  } else {
    pp$region_ref_position = ifelse(strand(pp)=="+", start(pp), end(pp) )
  }
  
  ovlp = findOverlaps(query_region, pp, maxgap=maxgap, ignore.strand=T)
  
  subject_hits_tmp = pp[subjectHits(ovlp)]
  
  # for palindromic target motifs we ignore any hits on the minus strand
  palindromic_not_on_ref_strand = (strand(subject_hits_tmp)!=palindrome_ref_strand & subject_hits_tmp$palindromic)
  ovlp_orig = ovlp
  ovlp = ovlp[!palindromic_not_on_ref_strand]
  
  qH = queryHits(ovlp)
  
  selected_target_regions = pp[subjectHits(ovlp)]
  selected_query_pos = query_region[qH]
  
  motif_strand_dir = ifelse(strand(selected_target_regions)=="+", 1,-1)
  bp_excl_rel = motif_strand_dir *  (start(selected_query_pos) - selected_target_regions$region_ref_position)
  
  if (break_between_nucleotides) {
    increment = motif_strand_dir
    bp_incl_rel = ifelse(strand(selected_query_pos)=="+", bp_excl_rel+increment, bp_excl_rel-increment)
    
    bp_dist = (bp_excl_rel + bp_incl_rel)/2
  } else {
    bp_dist = bp_excl_rel
  }
  
  min_bp_dist_per_query = tapply(abs(bp_dist), qH, min)
  
  result = data.frame(query_index = qH, distance = bp_dist, subject_index= subjectHits(ovlp), stringsAsFactors = F)
  result$min_bp_dist = min_bp_dist_per_query[as.character(qH)]
  result = subset(result, min_bp_dist == abs(distance))
  return(result)
}


if(FALSE) {
  # Test data for distance assessment
  pattern_dict = DNAStringSet(c("AAATTT","AAAATT","AATATT"))
  palindromic = (pattern_dict==reverseComplement(pattern_dict))
  
  g_seq = DNAString("CGAAAAAATTAAGCAAATTTCCCCTTTAAAGGGGAATTTTGGG")
  
  I2Rp <- function(x) {ff = unlist(lapply(pattern_pos_plus, length)); g = GRanges(rep("Z", sum(ff)), unlist(x[ff>0]), ind=(1:length(x))[ff>0], strand=rep("+", sum(ff))); return(g) }
  I2Rn <- function(x) {ff = unlist(lapply(pattern_pos_plus, length)); ii = unlist(x[ff>0]); s = length(g_seq)-start(ii)+1; e = length(g_seq)-end(ii)+1;  g = GRanges(rep("Z", sum(ff)), IRanges(e,s), ind=(1:length(x))[ff>0], strand=rep("-", sum(ff))); return(g) }

  pattern_pos_plus = matchPDict(pattern_dict, g_seq)
  pattern_pos_neg = matchPDict(pattern_dict, reverseComplement(g_seq))
  
  pattern_pos_test = c(I2Rp(pattern_pos_plus), I2Rn(pattern_pos_neg))
  pattern_pos_test$pattern = as.character(pattern_dict)[as.integer(pattern_pos_test$ind)]
  pattern_pos_test$palindromic = palindromic[as.integer(pattern_pos_test$ind)]
  
  dsb_ranges = GRanges(rep("Z",9), IRanges(c(4,6,8,10,15,20,37,39,30), c(4,6,8,10,15,20,37,39,30)), strand=c("+","+","-","-","+","+","-","-","+"), sampleID=rep("FOO", 9) )
  
  rr = compute_distances(dsb_ranges, pattern_pos_test, break_between_nucleotides = T, palindrome_ref_strand = "-")
  
  tmp_dsb_ends = dsb_ranges
  tmp_dsb_ends$min_bp_dist = NA_real_
  tmp_dsb_ends[rr$query_index]$min_bp_dist = rr$distance
  tmp_dsb_ends$motif_strand = NA_character_
  tmp_dsb_ends[rr$query_index]$motif_strand = strand(pattern_pos_test[rr$subject_index])
  tmp_dsb_ends = tmp_dsb_ends[!is.na(tmp_dsb_ends$min_bp_dist)]
  bp_regions_test = tmp_dsb_ends
  
  tmp2 = data.table(sampleID=tmp_dsb_ends$sampleID, strand=as.character(strand(tmp_dsb_ends)), dist = tmp_dsb_ends$min_bp_dist)
  tmp3 = tmp2[, .(count=.N), by=c("sampleID","strand","dist")]
}
