library(data.table)
setwd("/data_genome1/public_data/Ganniakis_colon_cancer/")
data = fread("Giannakis_som_mutations.txt.csv", sep="\t", skip = 3)
ed = fread("Giannakis_sample_description.txt.csv", sep="\t", skip=4)

save(data, ed, file="Giannakis_colon_cancer_mutations.Rdata")
