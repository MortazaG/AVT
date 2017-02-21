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
parser.add_argument("-f", "--fasta", help="FASTA - Overview")
parser.add_argument("-b", "--bam", help="BAM - Overview")
parser.add_argument("-bcr", help="BAM Graph - Coverage per reference ")
parser.add_argument("-bcp", help="BAM Graph - Coverage per position")

# Parse the above arguments, so that they can be used in the script.
args = parser.parse_args()

def fetch_fasta(ff):

    '''
    Fetch and print necessary information from a FASTA file.
    FASTA filehandle is received through the ff argument.

    Receives FASTA filename through ff argument.
    Outputs information to txt file in results/
    '''

    # Open FASTA file
    with open('results/' + ff + '.txt', 'w') as save_fasta:

        # Use SeqIO from BipPython to parse fasta file.
        for fasta_record in SeqIO.parse(ff, 'fasta'):
            gc_content = SeqUtils.GC(fasta_record.seq)

            # Write information from FASTA file to txt file.
            save_fasta.write('ID: %s\r\n' % fasta_record.id)
            save_fasta.write('GC content: %.2f%%\r\n' % gc_content)

def bam_view(bf):

    '''
    Write necessary information from a BAM file to txt file.
    BAM filehandle is received through the bf argument.

    Receives pysams samfile as argument through bf.
    Outputs information to txt file in results/ folder.
    '''

    # Calculate the total length of the reference sequence.
    # Achieved by summing up the lengths of individual references.
    ref_size = [length for length in bf.lengths]
    ref_size = sum(ref_size)

    # Fetch amount of reads that are available.
    tot_reads = bf.count()
    map_reads = bf.mapped
    unmap_reads = bf.unmapped

    # Calculate the sum of the reads lengths
    sum_reads_lengths = [len(read.seq) for read in bf.fetch()]
    sum_reads_lengths = sum(sum_reads_lengths)

    # Define average read length
    av_read_length = sum_reads_lengths / float(map_reads)

    # Calculate average coverage for the reference genome,
    # using the average read length.
    av_coverage = (av_read_length * map_reads) / float(ref_size)

    # Write extracted information from the BAM file to txt file
    with open('results/' + args.bam + '.txt', 'w') as save_bam:

        save_bam.write('Reference size: %d\r\n' % ref_size)
        save_bam.write('Total number of reads: %d\r\n' % tot_reads)
        save_bam.write('Mapped reads: %d\r\n' % map_reads)
        save_bam.write('Unmapped reads: %d\r\n' % unmap_reads)
        save_bam.write('Average read length: %.2f\r\n' % av_read_length)
        save_bam.write('Average coverage: %.4f\r\n' % av_coverage)

def bam_cov_ref(bf):

    '''
    Produce 'coverage per reference' graph, following:
    coverage = (read_length*read_count)/reference_length
    We will be using the average read length for each reference

    Receives pysams samfile as argument through bf.
    Outputs graphs as pdf files in results/graphs/ folder.
    '''

    # Create lists to hold reference names and lengths
    refs = bf.references
    refs_lengths = bf.lengths

    # Initialize list which will contain coverage for each reference
    refs_coverage = []

    # Loop through references by index nr and calculate coverage
    for i in range(len(refs)):

        # Define sum of read length for reference, so that we can calculate
        # average read length. This is done using pysams fetch()
        sum_read_length = [len(read.seq) for read in bf.fetch(refs[i])]
        sum_read_length = sum(sum_read_length)

        # Create list to hold total reads for reference
        tot_reads = bf.count(refs[i])
        av_read_length = sum_read_length / tot_reads

        # Define the length of current reference
        ref_length = refs_lengths[i]

        # Calculate coverage for reference according to given formula
        ref_coverage = (tot_reads * av_read_length) / float(ref_length)
        refs_coverage.append(ref_coverage)

    # Produce numerical values for x axis
    refs_range = [r+1 for r in range(len(refs))]

    # Plot graph using matplotlib.pyplot
    plt.figure('cov_ref')
    plt.plot(refs_range, refs_coverage)

    plt.title('Coverage per reference\n')
    plt.ylabel('Mean Coverage')
    plt.xlabel('Reference')
    plt.xticks(refs_range, refs, rotation='vertical')
    plt.grid()

    plt.savefig('results/graphs/' + args.bcr + '_cov_ref.pdf', bbox_inches='tight')
    plt.close()

def bam_cov_pos(bf):

    '''
    Produce 'coverage per position' graph for each individual reference

    Receive pysams samfile as argument through bf.
    Outputs graph as a pdf file in results/graphs/ folder.
    '''

    # Create lists to hold reference names and lengths
    refs = bf.references

    # Loop through individual references
    for ref in refs:

        # Use pysam pileup to get coverage for each position in reference
        ref_pos = [p.pos for p in bf.pileup(ref)]
        ref_pos_coverage = [p.n for p in bf.pileup(ref)]

        # Plot graph using matplotlib
        plt.figure()
        plt.plot(ref_pos, ref_pos_coverage)

        plt.title('Coverage across reference ' + ref + '\n')
        plt.ylabel('Mean Coverage')
        plt.xlabel('Position (bp)')
        plt.grid()

        plt.savefig('results/graphs/' + args.bcp + ref + '_cov_pos.pdf', bbox_inches='tight')
        plt.close()

def fetch_bam(bf):

    '''
    BAM filename is received through the bf argument and is then
    checked for Index file. With the Index file, samfile is created
    as a variable and passed on as an argument to the respective funtion.
    '''

    # Open BAM file using pysam
    samfile = pysam.AlignmentFile(bf, 'rb')

    # Check if BAM file has an index file. If False, give option
    # to create index file.
    if samfile.has_index() == False:
        print '\nNo BAM index file could be found for %s' % (bf)

        bam_bai = raw_input('Produce index file now? [y/n]: ').lower()

        if bam_bai == 'y':
            pysam.index(bf)

        elif bam_bai == 'n':
            print '\nThis program can\'t be run without a BAM index file.\n'

        else:
            print '\nSomething wen\'t wrong...\n'

        exit()

    # Call the respective function with samfile as argument.
    if args.bam:
        bam_view(samfile)

    if args.bcr:
        bam_cov_ref(samfile)

    if args.bcp:
        bam_cov_pos(samfile)

    samfile.close()

def main():

    '''
    Main function reads in FASTA and BAM files by calling fetch_fasta()
    and fetch_bam() functions with the user input as argument.
    '''

    # Check if args.fasta is present and then call fetch_fasta with
    # args.fasta as argument.
    if args.fasta:
        fetch_fasta(args.fasta)

    # Check which of the bam arguments is present and call fetch_bam on it.
    if args.bam:
        fetch_bam(args.bam)

    if args.bcr:
        fetch_bam(args.bcr)

    if args.bcp:
        fetch_bam(args.bcp)

    # Future capability to be added
    if args.fasta and args.bam:
        print 'Sorry, can\'t do that yet...'

# Makes sure main() is only run when this script is called from
# from itself.
if __name__ == "__main__":

    import os

    # Create necessary folders if they don't already exist
    try:
        os.makedirs('results/graphs')

    except OSError:
        pass

    main()
