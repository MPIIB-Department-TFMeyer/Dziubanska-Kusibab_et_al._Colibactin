library(GenomicFeatures)
# define genomic bins

# Chrom info from human_g1k_v37_decoy.fasta.fai 
# chrominfo_tmp = read.table(gzfile("/data_genome2/References/Human/Sequences/Genome/human_g1k_v37/human_g1k_v37_decoy.fasta.fai.gz"), sep="\t", header=F, stringsAsFactors = F)
# chrominfo_tmp$chr = sapply(strsplit(chrominfo_tmp$V1, " "), function(x) x[1])
# 
# chrom_info = data.frame(stringsAsFactors = F, V1=chrominfo_tmp$chr, V2=chrominfo_tmp$V2)

chrom_info = read.table("./ChromInfo.txt", sep="\t", header=T)

basic_gtf_file = "./gencode.v28lift37.basic.annotation_fixed_chroms.gtf"
chrom_info_txdb = chrom_info
colnames(chrom_info_txdb)=c("chrom","length")
chrom_info_txdb$is_circular=ifelse(chrom_info_txdb$chrom=="MT",T,F)
Gencode_basicTxDB = makeTxDbFromGFF(basic_gtf_file, format="gtf", chrominfo=chrom_info_txdb, dataSource="GENCODEv28lift37 basic hg19", organism="Homo sapiens")
saveDb(Gencode_basicTxDB, file="GENCODEv28lift37_basic_hg19_TxDB.db")
