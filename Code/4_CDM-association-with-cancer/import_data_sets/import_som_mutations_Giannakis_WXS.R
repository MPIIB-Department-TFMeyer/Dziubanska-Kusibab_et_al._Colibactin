library(data.table)
#setwd("/path/to/this/folder")
data = fread("../../../Data/Giannakis_WXS_CRC/Giannakis_som_mutations.txt.csv", sep="\t", skip = 3)
ed = fread("../../../Data/Giannakis_WXS_CRC/Giannakis_sample_description.txt.csv", sep="\t", skip=4)

save(data, ed, file="../../../Data/Giannakis_WXS_CRC/Giannakis_colon_cancer_mutations.Rdata")
