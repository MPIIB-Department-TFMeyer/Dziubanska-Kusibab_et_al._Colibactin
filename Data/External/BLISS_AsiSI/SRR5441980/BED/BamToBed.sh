bamToBed -i ../BWA/AsiSI_induced_sorted.bam | awk | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,".",$5,$6;}' | gzip -c -9 > AsiSI_induced_positions.bed.gz
