# Script5_Filter3_excessive_individuals.R

## R Script to remove the SNPs covered in an excessive number of individuals

rm(list = ls()) #clears environment
options(stringsAsFactors = F, mc.cores = 2)

# Import the .mafs file.
data<-read.table("Good_coverage_SNPs.txt",header=T)

# Create a pdf - the size of the graph can be adjusted with width and height.
pdf("histogram_Nb_ind.pdf",width=6,height=4)

# Set the margin of the graph.
par(mar = par("mar") + c(0, 1, 0, 0))

# Plot the histogram.
hist(data$nInd, breaks=40, xlab= "Number of individuals",ylab="", main="",
	asp=NA, las=1,yaxt='n')

# Rename the y axis for more clarity - values need to be changed according to the data.
axis(2, at=c(0, 20000,40000,60000,80000,100000,120000),
	labels=c("0", "20,000","40,000","60,000","80,000","100,000","120,000"),
	col.axis="black", las=2)

# Add the y axis label.
mtext("Number of SNPs", side = 2, line = 4)
dev.off()
# The cut-off is defined by the tail of the distribution.
# Then the SNPs for which there is coverage for >n individuals are filtered out from the .mafs file. As an example we use a cut-off of 10 individuals here.

# Select the SNPs covered in 10 individuals or less.
data2<-data[which(data$nInd<=10),]

# Create a new file for the filtered data.
write.table(data2,file="SNP_10ind.txt",dec=".",sep="\t", row.names=F, col.names=T, quote=F) 
