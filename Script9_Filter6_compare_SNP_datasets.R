# Script9_Filter6_compare_SNP_datasets.R

## R-script to compare 2 .mafs files

rm(list = ls()) #clears environment
options(stringsAsFactors = F, mc.cores = 2)

# Import the two .mafs files.
data1<-read.table("Ppho_SNPs_max20_300bp_new.geno.txt",header=F) 
data2<-read.table("Ppho_SNPs_max20_300bp.geno.txt",header=F)
colnames(data1)<-c("chromo","position","major","minor","unknownEM","pu.EM","nInd")
colnames(data2)<-c("chromo","position","major","minor","unknownEM","pu.EM","nInd")

# To find the shared SNPs a single column containing both the scaffold (or contig or chromosome) and position numbers is created.

# For data1, create a new column with the scaffold ID and the position on the scaffold.
position1<-paste(data1$chromo,data1$position)

# Add this new column to the table.
data1<-cbind(position1,data1)

# Rename the columns.
colnames(data1)<-c("pos","chromo","position","major","minor","unknownEM","pu.EM","nInd")

# The same steps are undertaken for data2.
position2<-paste(data2$chromo,data2$position)
data2<-cbind(position2,data2)
colnames(data2)<-c("pos","chromo","position","major","minor","unknownEM","pu.EM","nInd")

# Find the SNPs that are found in both data files.
duplicates<- data1[is.element(data1[,1], data2[,1]), ]

# Count the number of shared SNPs.
Nb_shared_SNPs<-nrow(duplicates)

# Remove the first column that is not needed any more, rename the columns and write a new text file containing the shared SNPs by the two datasets. 
duplicates2<-duplicates[,2:8]
colnames(duplicates2)<-c("chromo","position","major","minor","unknownEM","pu.EM","nInd")
write.table(duplicates2,file="shared_SNPs.txt",dec=".",sep="\t", row.names=F, col.names=T, quote=F) 
