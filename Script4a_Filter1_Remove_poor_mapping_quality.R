# Script4a_Filter1_Remove_poor_mapping_quality.R

### R Script to remove the SNPs that are in a region of poor mapping quality
## 1. The list of SNPs (.mafs file) and the regions of poor mapping quality are imported in R.

rm(list = ls()) # clears environment
options(stringsAsFactors = F, mc.cores = 2)

# Import the list of SNPs.
snp <- read.table("genolike_data.mafs", header=TRUE)

# Import the list of sites with poor mapping quality.
poor <- read.table("poor_mapping_quality.txt", header=FALSE) 

## 2. In `poor`, remove all scaffolds (or contigs or chromosomes depending 
#     on the assembly level of the reference genome) that do not contain any SNP.

# First, get the list of scaffolds containing SNPs.
scaffolds.snp <- unique(snp[,1])

# Then, get the list of scaffolds containing poor mapping quality regions.
scaffolds.poor <- unique(poor[,1])

# Identify which scaffolds with poor mapping quality do not contain any SNPs.
addr <- !is.element(scaffolds.poor, scaffolds.snp)
missing.scaffolds <- scaffolds.poor[addr]

# Remove the corresponding lines from `poor`: roughly one third of `poor` lines are 
# removed. This should speed up the process.
poor.clean <- poor[!is.element(poor[,1], missing.scaffolds), ]
poor.clean <- droplevels(poor.clean)
poor.clean <- poor.clean[order(poor.clean[,1]),]

## 3. In `snp`, remove all scaffolds that do not appear in `poor.clean`, 
#     i.e. scaffolds that contain SNPs and do not show "poor mapping quality" regions.
 
# Identify which scaffolds contain SNPs and are not present in `poor.clean`.
addr2 <- !is.element(unique(snp[,1]), unique(poor.clean[,1]))
snp.scaffolds.not.in.poor <- unique(snp[,1])[addr2]

# Remove the corresponding lines from snp.
snp.clean <- snp[ !is.element(snp[,1], snp.scaffolds.not.in.poor), ]
snp.clean <- droplevels(snp.clean)
snp.clean <- snp.clean[order(snp.clean[,1]),]

# Store the lines removed in a new object called snp.ok. These SNPs will be added 
# later to the list of SNPs that are kept.
snp.ok <- snp[ is.element(snp[,1], snp.scaffolds.not.in.poor), ]

# We now have 2 objects containing the same list of scaffolds, sorted in the same order:
identical(unique(snp.clean[,1]), unique(poor.clean[,1]))
# The scaffolds in those objects contain SNPs and at least some regions have a "poor mapping quality".

##4.  The SNPs from `snp.clean` that fall within a poor mapping quality region need to be filtered out.
# Since `snp.clean` and `poor.clean` are smaller than the raw datasets `snp` and `poor` respectively,
# the following computations should be slightly less time-consuming than with the original datasets.
nrow(snp) ; nrow(snp.clean)
nrow(poor) ; nrow(poor.clean)

# First, get the list of SNPs for each scaffold in snp.clean.
snp.list <- tapply(snp.clean$position, snp.clean$chromo, as.vector)

# Then, for each scaffold, get a list of regions for which there is poor mapping quality.
library(plyr) # library plyr needs to be installed
tmp1 <- dlply(poor.clean[,2:3], .(poor.clean[,1]))
poor.list <- lapply(tmp1, function(x) as.vector(unlist(apply(x,1, function(y) y[1]:y[2]))))

# Check that the scaffolds in both lists are the same and that they are both ordered the same way.
if ( !identical(names(snp.list), names(poor.list)) ) {

    stop("The lists of scaffolds are not compatible")    
    
} else {
    # If so, create a function that returns the list of SNPs to keep for one scaffold.
    keep.snp <- function(stest, ptest) {
        !is.element(stest,ptest)
    }
    
    # Then apply this function to the complete list.
    # Each element of this list is a vector of TRUE or FALSE: TRUE when the SNP should be kept, 
    # FALSE otherwise.
    snp.to.keep <- mapply(keep.snp, snp.list, poor.list)

    # Convert the format from list to vector.
    names(snp.to.keep) <- NULL
    snp.to.keep <- unlist(snp.to.keep)
    
    # This logical vector is used to create two new tables containing the SNPs that should be 
    # kept and the SNPs that fall within a poor mapping quality region.
    kept.snp <- snp.clean[snp.to.keep, ]
    dropped.snp <- snp.clean[!snp.to.keep, ]
}

## 5. Finally, the object `snp.ok` that contains the SNPs in scaffolds that have no region of 
#     poor mapping quality is added to `kept.snp`.
kept.snp <- rbind(kept.snp, snp.ok)
kept.snp <- kept.snp[order(kept.snp[,1]),]

# A file is written for the selected SNPs
write.table(kept.snp, file = "SNPs_goodquality.txt", row.names = FALSE, sep = "\t",quote=F)
