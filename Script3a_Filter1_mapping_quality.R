## Rscript3a_Filter1_mapping_quality.R
## R Script to select the poor mapping quality regions ("Filter 1")

rm(list = ls()) #clears environment
options(stringsAsFactors = F, mc.cores = 2)

# Import the data.
data<-read.table("data_callable.bed", header=F)

# Select the regions of poor mapping quality.
data2<-data[which(data[,4]=="POOR_MAPPING_QUALITY"),] 

# Produce a new file containing only the regions of poor mapping quality.
write.table(data2,file="poor_mapping_quality.txt",dec=".",sep="\t", row.names=F, col.names=F, quote=F) 
