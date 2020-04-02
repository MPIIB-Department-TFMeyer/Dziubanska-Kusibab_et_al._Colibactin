GENOME_INDEX=/data_genome1/References/GenomeIndices/combined/Human_CtrL2_Ecoli/BWA/human_g1k_v37_CtrL2_EcoliM1.fasta
bwa mem -t 8 $GENOME_INDEX ../FASTQ/AsiSI_induced_genomic.fastq.gz | samtools view -bu | samtools sort -@ 4 -T tmp -O BAM -o AsiSI_induced_sorted.bam
