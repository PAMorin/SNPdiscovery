# Script4b_Filter2_Remove_excessive_coverage.R

### R Script to remove the SNPs that are in a region of excessive coverage
## 1. The list of SNPs (output file of Filter 1 - script 4a) and the regions of excessive coverage are imported in R.

rm(list = ls()) #clears environment
options(stringsAsFactors = F, mc.cores = 2)

# Import the list of SNPs.
snp <- read.table("SNPs_goodquality.txt", header=TRUE)

# Import the list of sites with excessive coverage.
high_cov <- read.table("excessive_coverage.txt", header=FALSE) 

## 2. In `high_cov`, remove all scaffolds (or contigs or chromosomes depending
#     on the assembly level of the reference genome) that do not contain any SNP.

# First, get the list of scaffolds containing SNPs.
scaffolds.snp <- unique(snp[,1])

# Then, get the list of scaffolds containing excessive coverage regions.
scaffolds.high_cov <- unique(high_cov[,1])

# Identify which scaffolds with excessive coverage do not contain any SNPs.
addr <- !is.element(scaffolds.high_cov, scaffolds.snp)
missing.scaffolds <- scaffolds.high_cov[addr]

# Remove the corresponding lines from `high_cov`: roughly one third of `high_cov` lines are
# removed. This should speed up the process.
high_cov.clean <- high_cov[!is.element(high_cov[,1], missing.scaffolds), ]
high_cov.clean <- droplevels(high_cov.clean)
high_cov.clean <- high_cov.clean[order(high_cov.clean[,1]),]

## 3. In `snp`, remove all scaffolds that do not appear in `high_cov.clean`, 
#     i.e. scaffolds that contain SNPs and do not show "excessive coverage" regions.
 
# Identify which scaffolds contain SNPs and are not present in `high_cov.clean`.
addr2 <- !is.element(unique(snp[,1]), unique(high_cov.clean[,1]))
snp.scaffolds.not.in.high_cov <- unique(snp[,1])[addr2]

# Remove the corresponding lines from snp.
snp.clean <- snp[ !is.element(snp[,1], snp.scaffolds.not.in.high_cov), ]
snp.clean <- droplevels(snp.clean)
snp.clean <- snp.clean[order(snp.clean[,1]),]

# Store the lines removed in a new object called snp.ok. These SNPs will be added
# later to the list of SNPs that are kept.
snp.ok <- snp[ is.element(snp[,1], snp.scaffolds.not.in.high_cov), ]

# We now have 2 objects containing the same list of scaffolds, sorted in the same order:
identical(unique(snp.clean[,1]), unique(high_cov.clean[,1]))
# The scaffolds in those objects contain SNPs and at least some regions have an "excessive coverage".

##4.  The SNPs from `snp.clean` that fall within an excessive coverage region need to be filtered out.
# Since `snp.clean` and `high_cov.clean` are smaller than the raw datasets `snp` and `high_cov` respectively, 
# the following computations should be slightly less time-consuming than with the original datasets.
nrow(snp) ; nrow(snp.clean)
nrow(high_cov) ; nrow(high_cov.clean)

# First, get the list of SNPs for each scaffold in snp.clean.
snp.list <- tapply(snp.clean$position, snp.clean$chromo, as.vector)

# Then, for each scaffold, get a list of regions for which there is excessive coverage.
library(plyr) # library plyr needs to be installed
tmp1 <- dlply(high_cov.clean[,2:3], .(high_cov.clean[,1]))
high_cov.list <- lapply(tmp1, function(x) as.vector(unlist(apply(x,1, function(y) y[1]:y[2]))))

# Check that the scaffolds in both lists are the same and that they are both ordered the same way.
if ( !identical(names(snp.list), names(high_cov.list)) ) {

    stop("The lists of scaffolds are not compatible")    
    
} else {
    # If so, create a function that returns the list of SNPs to keep for one scaffold.
    keep.snp <- function(stest, ptest) {
        !is.element(stest,ptest)
    }
    
    # Then apply this function to the complete list.
    # Each element of this list is a vector of TRUE or FALSE: TRUE when the SNP should be kept, 
    # FALSE otherwise.
    snp.to.keep <- mapply(keep.snp, snp.list, high_cov.list)

    # Convert the format from list to vector.
    names(snp.to.keep) <- NULL
    snp.to.keep <- unlist(snp.to.keep)
    
    # This logical vector is used to create two new tables containing the SNPs that should be 
    # kept and the SNPs that fall within an excessive coverage region.
    kept.snp <- snp.clean[snp.to.keep, ]
    dropped.snp <- snp.clean[!snp.to.keep, ]
}

## 5. Finally, the object `snp.ok` that contains the SNPs in scaffolds that have no region of 
#     excessive coverage is added to `kept.snp`.
kept.snp <- rbind(kept.snp, snp.ok)
kept.snp <- kept.snp[order(kept.snp[,1]),]

# A file is written for both objects.
write.table(kept.snp, file = "Good_coverage_SNPs.txt", row.names = FALSE, sep = "\t",quote=F)
