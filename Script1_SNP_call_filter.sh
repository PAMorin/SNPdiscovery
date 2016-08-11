####
# Shell script to perform samtools pileup, call SNPs with bcftools, and filter SNPs for depth, quality, and flanking region.
# Requires several programs (see below) and calls python script "Script2_generate_genotype_blocks.py" that must also be in same directory.
#Call script with format ./script.sh file.bam file.fa 300 output
# the fragment size "300" is the default value for the total fragment size to be screened
#   (150bp on either side of the target SNP). This can be changed if larger or smaller
#   flanking regions are required.
# bam file must be sorted
# Reference fasta file must be indexed for samtools

input=`basename $1 .bam`
reference=$2
fragsize=$3
output=$4

# Perform a pileup 
# uses samtools; provide path/samtools if not in the system environment path.
samtools mpileup -g -f $reference $input.bam > $input_raw.bcf

# Call SNPS, and filter out any depth of coverages below and above given values (-d, -D).
# Uses bcftools; provide path/bcftools if not in the system environment path.
bcftools call -c -v  $input_raw.bcf | vcfutils.pl varFilter -d 10 -D 150 > $input_flt.vcf
# varfilter -d 10 -D 100 filters SNPs based on minimum coverage of 10 and maximum of 100. 

# Run custom Python script to generate SNP blocks.
# Uses python; provide path/python if not in the system environment path.
python Script2_generate_genotype_blocks.py -d $fragsize -v $input_flt.vcf -f $reference -q 50 > $output.geno.fasta 2> $output.geno.vcf
# -d 300  specifies that no SNPs closer than 150bp on either side are identified, and it will print out 300/2=150bp on either side of the SNPs in the .fasta output file. 
####
