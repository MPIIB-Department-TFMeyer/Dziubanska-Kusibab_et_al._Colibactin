suppressPackageStartupMessages(library('rtracklayer'))
library(naturalsort)
library(GenomeInfoDb)

# download the following file from the following URL
#url="ftp://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign100mer.bw"
url="wgEncodeCrgMapabilityAlign100mer.bw"
bwf = BigWigFile(url)
mappability_track_100mer=import(bwf)
seqlevelsStyle(mappability_track_100mer) <- "NCBI"

###############################################################
#### Bins 
###############################################################
chrom_info = read.table("hg19.genome",sep="\t",header=F, stringsAsFactors = F)
colnames(chrom_info) = c("chrom","size")
chrom_info_fixed = chrom_info
chrom_info_fixed = chrom_info_fixed[naturalorder(chrom_info_fixed$chrom),]
seqinfo = Seqinfo(seqnames=chrom_info_fixed$chrom, seqlengths=chrom_info_fixed$size, isCircular=rep(F, nrow(chrom_info_fixed)), genome="hg19")

bin_size=300
bin_starts_per_chrom = tapply(chrom_info_fixed$size, factor(chrom_info_fixed$chrom, levels=chrom_info_fixed$chrom), function(x) seq(1,x,bin_size),simplify=F)
c_names = do.call(c, mapply(names(bin_starts_per_chrom), lapply(bin_starts_per_chrom, length), FUN=function(c,n) return(rep(c,n)) ))
c_starts = do.call(c, bin_starts_per_chrom)
c_ends = c_starts + bin_size-1
last_bin_per_chrom = cumsum(lapply(bin_starts_per_chrom, length))
c_ends[last_bin_per_chrom] = chrom_info_fixed$size
genomic_bins_300bp = GRanges(c_names, IRanges(c_starts, c_ends), seqinfo=seqinfo)
# # Make sure that the chromosomes in chrom_info have the same sort order as the ones in genomic_bins, otherwise the result might be invalid
###############################################################

genomic_bins_300bp$avg_mappability_score = 0
for (chr in seqnames(seqinfo(genomic_bins_300bp)) ) {
  c_bins = genomic_bins_300bp[seqnames(genomic_bins_300bp) == chr]
  bin_width = width(c_bins[1])
  c_mappability = mappability_track_100mer[seqnames(mappability_track_100mer) == chr]
  g = gaps(c_mappability)
  g_fixed = g[seqnames(g)==chr & strand(g)=="*"]
  g_fixed$score=0
  c_mappability_fixed = sort(c(c_mappability, g_fixed))
  mappability_expanded = Rle(rep(c_mappability_fixed$score, times=width(c_mappability_fixed)))
  # we here (ab)use the running mean to compute mean scores for (bin_width) sized intervals and then pich the center one
  bin_centers = seq(1,max(end(c_bins)), by=bin_width)+round(bin_width/2)-1
  bin_centers[length(bin_centers)] = max(end(c_bins))-1
  avg_per_bin = as.numeric(runmean(mappability_expanded, k=bin_width, endrule="constant")[bin_centers])
  genomic_bins_300bp[seqnames(genomic_bins_300bp) == chr]$avg_mappability_score = avg_per_bin
}

low_mappability_regions = reduce(genomic_bins_300bp[genomic_bins_300bp$avg_mappability_score < 0.9])
low_map_regions_df = data.frame(chrom=seqnames(low_mappability_regions), chrom_start=start(low_mappability_regions), chrom_end = end(low_mappability_regions))

#write.table(low_map_regions_df, file = "low_mappability_hg19__100mer_300bp_bins.bed",sep="\t", col.names=F, quote=F, row.names = F)

####################################################

save(low_mappability_regions, file="low_mappability_100mer_h19.Rdata")


