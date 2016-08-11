#! /bin/env python

import os, sys
from optparse import OptionParser
from itertools import groupby

def load_contigs_into_dict(options):
    """ Taken from: https://www.biostars.org/p/710/#1412"""

    # Load all the contig entries into a dictionary.
    contigs = {}

    with open(options.reference_filename, 'r') as reference:
        faiter = (x[1] for x in groupby(reference, lambda line: line[0] == ">"))
        for header in faiter:
            # drop the ">"
            header = header.next()[1:].strip().split()[0]
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())
            contigs[header] = seq

    return contigs


def load_vcf_entries_into_list(options):

    """ Example VCF file
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  bwa.bwa.sorted.bam
    scaffold1|size664435    2460    .       G       C       44.0073 .       DP=13;VDB=0.244943;SGB=-0.556411;RPB=0.988166;MQB=0.785785;MQSB=0.935229;BQB=0.97357;MQ0F=0;AF1=0.5;AC1=1;DP4=6,3,2,2;MQ=36;FQ=47.0153;PV4=1,0.364896,0.069705,1        GT:PL   0/1:74,0,180
    scaffold1|size664435    2669    .       TGG     TG      31.4691 .       INDEL;IDV=4;IMF=0.25;DP=16;VDB=0.414899;SGB=-0.556411;MQSB=1;MQ0F=0;AF1=0.5;AC1=1;DP4=10,1,0,4;MQ=37;FQ=34.4775;PV4=0.003663,1.87112e-11,1,1    GT:PL   0/1:69,0,182
    scaffold1|size664435    4715    .       A       G       15.2078 .       DP=3;VDB=0.06;SGB=-0.453602;RPB=1;MQB=1;BQB=1;MQ0F=0;AF1=0.508025;AC1=1;DP4=0,1,0,2;MQ=37;FQ=-12.2562;PV4=1,1,1,1       GT:PL   0/1:45,0,15
    scaffold1|size664435    4850    .       T       C       25.0206 .       DP=7;VDB=0.0707707;SGB=-0.556411;RPB=0.607698;MQB=1.01283;MQSB=1;BQB=1.01283;MQ0F=0;AF1=0.500002;AC1=1;DP4=1,2,0,4;MQ=37;FQ=21.0757;PV4=0.428571,0.340975,1,1   GT:PL   0/1:55,0,49
    scaffold1|size664435    5201    .       T       C       95.0077 .       DP=17;VDB=0.629475;SGB=-0.636426;RPB=0.447194;MQB=1;MQSB=1;BQB=0.447194;MQ0F=0;AF1=0.5;AC1=1;DP4=4,6,5,2;MQ=37;FQ=98.0159;PV4=0.334842,0.0991446,1,1    GT:PL   0/1:125,0,192
    scaffold1|size664435    5392    .       C       G       86.0076 .       DP=12;VDB=0.0328229;SGB=-0.616816;RPB=0.304946;MQB=0.864013;MQSB=0.896474;BQB=0.71154;MQ0F=0;AF1=0.5;AC1=1;DP4=3,3,4,2;MQ=36;FQ=82.8224;PV4=1,1,1,1     GT:PL   0/1:116,0,111
    scaffold1|size664435    5398    .       A       G       97.0078 .       DP=12;VDB=0.031576;SGB=-0.616816;RPB=0.304946;MQB=0.864013;MQSB=0.896474;BQB=0.71154;MQ0F=0;AF1=0.5;AC1=1;DP4=3,3,4,2;MQ=36;FQ=96.4769;PV4=1,0.335516,1,1       GT:PL   0/1:127,0,126
    scaffold1|size664435    5524    .       CTT     CT      32.4719 .       INDEL;IDV=3;IMF=0.333333;DP=9;VDB=0.200442;SGB=-0.511536;MQSB=1;MQ0F=0;AF1=0.500396;AC1=1;DP4=0,2,3,0;MQ=37;FQ=-7.39771;PV4=0.1,1,1,0.345076    GT:PL   0/1:70,0,28
    """

    entries = []
    with open(options.vcf_filename, 'r') as vcf_file:

        for line in vcf_file:
            if line.startswith("#"):
                continue
            else:
                entries.append(line.strip().split())

    return entries


def write_genotypes(options, contigs, ref_contig, position, reference, alternate):

    sequence = contigs[ref_contig]
    position = int(position)

    if 'N' in sequence[position-(options.distance/2)-1:position+(options.distance/2)] or \
            len(sequence[position-(options.distance/2)-1:position+(options.distance/2)]) < options.distance + 1:
        return False

    print '>' + ref_contig + ',pos=' + str(position) + ',geno=' + reference + '/' + alternate
    print sequence[position-(options.distance/2)-1:position+(options.distance/2)]

#    print '>' + ref_contig + ',pos=' + str(position) + ',geno=' + alternate
#    print sequence[position-(options.distance/2)-1:position-1] + alternate + sequence[position:position+(options.distance/2)]

    return True


def get_options():
    parser = OptionParser()

    # Input parameters.
    parser.add_option("-v", "--vcf", dest="vcf_filename", help="VCF filename.")
    parser.add_option("-f", "--reference", dest="reference_filename", help="Reference FASTA filename.")
    parser.add_option("-q", "--quality", dest="quality", default=50.0, type=float)
    parser.add_option("-d", "--distance", dest="distance", default=300, type=int)

    (options, args) = parser.parse_args()

    return (options,args)


def main():

    (options, args) = get_options()

    # Load all the contig entries into a dictionary.
    contigs = load_contigs_into_dict(options)

    # Load all the VCF entries into a dictionary.
    entries = load_vcf_entries_into_list(options)

    for i in range(0, len(entries)):

        # Does the current SNP have a quality greater than or equal to the required value?
        if float(entries[i][5]) >= options.quality and 'INDEL' not in entries[i][7]:

            # Are we within an acceptable distance from the previous SNP?
            if i == 0 or entries[i-1][0] != entries[i][0] or \
                    int(entries[i][1]) - int(entries[i-1][1]) > options.distance:

                # Are we within an acceptable distance of the next SNP?
                if i == len(entries)-1 or entries[i+1][0] != entries[i][0] or \
                        int(entries[i+1][1]) - int(entries[i][1]) > options.distance:
                    if write_genotypes(options, contigs, entries[i][0], entries[i][1], entries[i][3], entries[i][4]):
                                    sys.stderr.write('\t'.join(entries[i]) + '\n')


if __name__ == '__main__':
    main()
