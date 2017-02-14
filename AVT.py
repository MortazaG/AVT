#!/usr/bin/env python

# Licence:

# Usage:

import sys
import pysam
import matplotlib.pyplot as plt
from Bio import SeqIO, SeqUtils

# Check if user has argparse module available.
# Oterwise, raise an error and ask the user if it should be installed.
try:
    import argparse
except ImportError:
	sys.stderr.write("[Error] The python module 'argparse' is not installed\n")
	sys.stderr.write("[--] Would you like to install it now using 'sudo easy_install' [Y/N]? ")
	answer = sys.stdin.readline()
	if answer[0].lower() == "y":
		sys.stderr.write("[--] Running 'sudo easy_install argparse'\n")
		from subprocess import call
		call(["sudo", "easy_install", "argparse"])
	else:
		sys.exit("[Error] Exiting due to missing dependency 'argparser'")

# Initialization of argparser, with the name of the program and
# the necessary arguments.
parser = argparse.ArgumentParser(prog="avt.py")
parser.add_argument("-f", "--fasta", help="Enter name of FASTA file")
parser.add_argument("-b", "--bam", help="Enter name of BAM file")

# Parse the above arguments, so that they can be used further down.
args = parser.parse_args()

def main():

    '''
    Main function reads in FASTA and BAM files in order to produce
    necessary information from them.
    '''

    # Check if argument 'args.fasta' is present and print the required
    # information from the FASTA file.
    if args.fasta:
        for fasta_record in SeqIO.parse(args.fasta, 'fasta'):
            fasta_id = fasta_record.id
            gc_content = SeqUtils.GC(fasta_record.seq)

            print 'ID: %s' % fasta_id
            print 'GC content: %.2f%%\n' % gc_content

    # Check if argument 'args.bam' is present and open it op for reading.
    if args.bam:
        samfile = pysam.AlignmentFile(args.bam, 'rb')

        # Check if BAM file has an index file. If False, give info about
        # how to create index file using samtools.
        # Alternatively we could use check_index() method as well.
        if samfile.has_index() == False:
            print '\nNo BAM index file could be found for %s.\
                    \nUse \'samtools index %s\' to create an index file. \
                    \nThe resulting file should be named as such: %s.bai\n'\
                    % (args.bam, args.bam, args.bam)

        else:

            # Calculate the total length of the reference sequence.
            # Achieved by summing up the lengths of individual references,
            # such as chromosomes.
            ref_size = 0
            for length in samfile.lengths:
                ref_size += length

            # Fetch amount of reads that are available.
            tot_reads = samfile.count()
            map_reads = samfile.mapped
            unmap_reads = samfile.unmapped

            # Calculate the sum of the read lengths and
            # the GC content of each read.
            sum_read_length = 0
            gc_read = []
            for read in samfile.fetch():
                sum_read_length += len(read.seq)
                gc_read.append(SeqUtils.GC(read.seq))

            # Define average read length
            av_read_length = sum_read_length / float(map_reads)

            # Calculate average coverage for the reference genome,
            # using the average read length.
            av_coverage = (av_read_length * map_reads) / float(ref_size)

            # Print out all the available information from the BAM file.
            print '\nReference size: %d' % ref_size
            print 'Total number of reads: %d' % tot_reads
            print 'Mapped reads: %d' % map_reads
            print 'Unmapped reads: %d' % unmap_reads
            print 'Average read length: %.2f\n' % av_read_length

            print 'Average coverage: %.4f\n' % av_coverage

        samfile.close()

if __name__ == "__main__":
    main()
