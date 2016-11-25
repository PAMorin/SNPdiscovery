# Script11_Filter7_Remove_HWEexcessHet.R

### R Script to remove the SNPs that are in a region of excessive coverage
## 1. The list of SNPs (output file of Filter 1 - script 4a) and the regions of excessive coverage are imported in R.

rm(list = ls()) #clears environment
options(stringsAsFactors = F, mc.cores = 2)

# Import the list of SNPs.
snp <- read.table("SNPs_MAF_good.txt", header=TRUE) #was header=TRUE, but that made row 2 the header

# Import the list of sites with excessive coverage.
high_het <- read.table("SNPs_hweExcessHet.txt", header=TRUE) 

##2.  The SNPs from `snp.clean` that fall within an excessive heterozygosity region need to be filtered out.
# Since `snp.clean` and `high_het.clean` are smaller than the raw datasets `snp` and `high_het` respectively, 
# the following computations should be slightly less time-consuming than with the original datasets.
nrow(snp) # ; nrow(snp.clean)
nrow(high_het) # ; nrow(high_het.clean)
# remove semicolon and nrow() if sections 2 and 3 deleted

# make vectors of concatenated chromo, position from each dataframe
snp.list <- paste(snp$chromo, "_", snp$position, sep="")
high_het.list <- paste(high_het$chromo, "_", high_het$position, sep="")

#determine which values are in both lists (vector of TRUE and FALSE, length of snp)
joint.list <- snp.list %in% high_het.list

#extract vector of rows that are not in common ("FALSE" in joint.list)
row_to_keep = which(joint.list == FALSE) 
#row_to_remove = which(joint.list == TRUE)

#new files including only the rows from snp that are not in high_het
kept.snp <- snp[row_to_keep,]

# write new file with kept SNPs.
write.table(kept.snp, file = "good_hwe_SNPs.txt", row.names = FALSE, sep = "\t",quote=F)
