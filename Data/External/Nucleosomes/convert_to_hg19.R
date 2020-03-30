library(data.table)
d = fread("zcat positioning_scores/peaks.min_peak_score_0.6.thresh_0.5.txt.gz", sep="\t", header=F)

# positions are 1-based according to readme.txt, convert to BED conventions (0-based)
d_new = data.table(chrom=d$V1, start=as.integer(d$V4-1), end=as.integer(d$V4))

fwrite(d_new, file="positioning_scores/peaks.min_peak_score_0.6.thresh_0.5_midpoints_for_liftover.txt", sep="\t", col.names = F)

# Download hg18 to hg19 chain files from UCSC
system("/data_genome1/References/Human/Sequences/LiftOver/liftOver positioning_scores/peaks.min_peak_score_0.6.thresh_0.5_midpoints_for_liftover.txt hg18ToHg19.over.chain.gz positioning_scores/peaks.min_peak_score_0.6.thresh_0.5_midpoints_hg19.txt positioning_scores/unmapped.txt")
