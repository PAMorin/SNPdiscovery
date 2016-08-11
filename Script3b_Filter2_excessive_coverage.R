## Rscript3b_Filter2_excessive_coverage
## R Script to select the excessive coverage regions ("Filter 2")

rm(list = ls()) #clears environment
options(stringsAsFactors = F, mc.cores = 2)

# Import the data.
data<-read.table("data_callable.bed", header=F)

# Select the regions of excessive coverage.
data3<-data[which(data[,4]=="EXCESSIVE_COVERAGE"),] 

# Produce a new file containing only the regions of excessive coverage.
write.table(data3,file="excessive_coverage.txt",dec=".",sep="\t", row.names=F, col.names=F, quote=F) 
