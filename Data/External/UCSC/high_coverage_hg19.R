# define SeqDups sites
library(GenomicRanges)

# input tables downloaded from UCSC table browser
HighCov_site_file="hg19_hiSeqDepthTop_5_Pct.txt.gz"
HighCov_site_tab = read.table(gzfile(HighCov_site_file), sep="\t",header=T,dec=".", comment.char="", as.is=T)
HighCov_site_tab$chrom = gsub("chr","",HighCov_site_tab$chrom)
HighCov_site_ranges_tmp = GRanges(HighCov_site_tab$chrom, IRanges(HighCov_site_tab$chromStart, HighCov_site_tab$chromEnd))
HighCov_site_ranges_5pct = reduce(HighCov_site_ranges_tmp)

HighCov_site_file="hg19_hiSeqDepthTop_1_Pct.txt.gz"
HighCov_site_tab = read.table(gzfile(HighCov_site_file), sep="\t",header=T,dec=".", comment.char="", as.is=T)
HighCov_site_tab$chrom = gsub("chr","",HighCov_site_tab$chrom)
HighCov_site_ranges_tmp = GRanges(HighCov_site_tab$chrom, IRanges(HighCov_site_tab$chromStart, HighCov_site_tab$chromEnd))
HighCov_site_ranges_1pct = reduce(HighCov_site_ranges_tmp)
save(HighCov_site_ranges_5pct, HighCov_site_ranges_1pct, file="HighCov_ranges.Rdata")
