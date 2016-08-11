#! /bin/env python

#SNP sequence selection from .maf file based on specified number of nucleotides surrounding the SNP (default = 300), and any N's within the flanking sequence. 

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

    """ Example mafs file
       chromo	position	major	minor	unknownEM	pu.EM	nInd
		NW_004197515.1	96746	T	A	0.249928	2.545278e-07	4
		NW_004197515.1	185333	T	C	0.181731	5.061095e-07	5
		NW_004197515.1	196007	T	C	0.499422	4.6215e-07	4
		NW_004197515.1	199416	G	A	0.181716	6.426417e-07	5
		NW_004197515.1	393950	T	C	0.182084	8.147291e-07	4
		NW_004197515.1	438598	T	C	0.499927	9.463028e-07	4
		NW_004197515.1	523279	A	G	0.499544	4.621587e-07	4
		NW_004197515.1	598251	C	T	0.499823	1.780898e-07	4
		NW_004197515.1	598857	T	C	0.249938	2.005634e-07	4
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
#    parser.add_option("-q", "--quality", dest="quality", default=50.0, type=float)
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
#        if float(entries[i][5]) >= options.quality and 'INDEL' not in entries[i][7]:

            # Are we within an acceptable distance from the previous SNP?
            if i == 0 or entries[i-1][0] != entries[i][0] or \
                    int(entries[i][1]) - int(entries[i-1][1]) > options.distance:

                # Are we within an acceptable distance of the next SNP?
                if i == len(entries)-1 or entries[i+1][0] != entries[i][0] or \
                        int(entries[i+1][1]) - int(entries[i][1]) > options.distance:
                    if write_genotypes(options, contigs, entries[i][0], entries[i][1], entries[i][2], entries[i][3]):
                                    sys.stderr.write('\t'.join(entries[i]) + '\n')


if __name__ == '__main__':
    main()
