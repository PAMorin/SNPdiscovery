####
# Shell script to extract flanking regions for SNPs, filtered to exclude sequences with N's
# Calls python script "Script2_generate_genotype_blocks.py" that must also be in same directory.
#Call script with format 
#./Script7_SNP_call_filter_GATK.sh Test_Ppho_SNPs_MAF_good.txt Ppho2Npho_consensus.fasta 300 output
# Reference fasta file must be indexed for samtools

input=`basename $1 .txt`
reference=$2
fragsize=$3
output=$4

# Run custom Python script to generate SNP blocks.
# Uses python; provide path/python if not in the system environment path.
python Script8_generate_genotype_blocks_GATK.py -d $fragsize -v $input.txt -f $reference > $output.geno.fasta 2> $output.geno.txt
# -d 300  specifies that no SNPs closer than 150bp on either side are identified, and it will print out 300/2=150bp on either side of the SNPs in the .fasta output file. 
####
