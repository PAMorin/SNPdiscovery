# Script6_Filter4_Remove_rare_SNPs.R

## R Script to select SNPs with minor allele frequency (MAF) >=0.05

rm(list = ls()) #clears environment
options(stringsAsFactors = F, mc.cores = 2)

# Download the data
data<-read.table("SNP_10ind.txt",header=T)

# select the SNPs with MAF >= 0.05
data2<-data[which(data$unknownEM>=0.05),]

# Create a file with the SNPs with MAF >= 0.05
colnames(data2)<-c("#chromo","position","major","minor","unknownEM","pu.EM","nInd")

write.table(data2,file="SNPs_MAF_good.txt",dec=".",sep="\t", row.names=F,col.names=T,quote=F) 
