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

# Parse the above arguments, so that they can be used in the script.
args = parser.parse_args()

def fetch_fasta(ff):

    '''
    Fetch and print necessary information from a FASTA file.
    FASTA filehandle is received through the ff argument.
    '''

    # Use SeqIO from BipPython to parse fasta file.
    for fasta_record in SeqIO.parse(ff, 'fasta'):
        gc_content = SeqUtils.GC(fasta_record.seq)

        print 'ID: %s' % fasta_record.id
        print 'GC content: %.2f%%\n' % gc_content

def fetch_bam(bf):

    '''
    Fetch and print necessary information from a BAM file.
    BAM filehandle is received through the bf argument.
    '''

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

        print 'Average coverage: %.4f' % av_coverage
        # print 'GC content: %.2f%%\n' % (sum(gc_read)/len(gc_read))

        bam_graphs(samfile)

        samfile.close()

def bam_graphs(sf):
    '''
    Produce necessary graphs for BAM file.
    - Coverage per reference
    - Coverage per position for each reference
    '''

    # Produce 'coverage per reference' graph
    refs = sf.references
    refs_lengths = sf.lengths

    refs_range = [r+1 for r in range(len(refs))]
    refs_coverage = []

    for i in range(len(refs)):

        sum_read_length = 0
        amt_reads = 0
        for read in sf.fetch(refs[i]):
            sum_read_length += len(read.seq)
            amt_reads += 1

        av_read_length = sum_read_length / amt_reads
        tot_reads = sf.count(refs[i])
        ref_length = refs_lengths[i]

        ref_coverage = (tot_reads * av_read_length) / float(ref_length)
        refs_coverage.append(ref_coverage)

    plt.figure('cov_ref')
    plt.plot(refs_range, refs_coverage)

    plt.title('Coverage per reference\n')
    plt.ylabel('Mean Coverage')
    plt.xlabel('Reference')
    plt.xticks(refs_range, refs, rotation='vertical')
    plt.grid()

    plt.savefig('graphs/' + args.bam + '.pdf', bbox_inches='tight')
    plt.close()

    # Produce 'coverage per position' graph for each reference
    # for i in range(len(refs)):
    for ref in refs:

        ref_pos = []
        ref_pos_coverage = []
        for b in sf.pileup(ref):
            ref_pos.append(b.pos)
            ref_pos_coverage.append(b.n)

        plt.figure(1)
        plt.plot(ref_pos, ref_pos_coverage)

        plt.title('Coverage across reference ' + ref + '\n')
        plt.ylabel('Mean Coverage')
        plt.xlabel('Position (bp)')
        plt.grid()

        plt.savefig('graphs/' + ref + '.pdf', bbox_inches='tight')
        plt.close()



def main():

    '''
    Main function reads in FASTA and BAM files by calling fetch_fasta()
    and fetch_bam() functions with the user input as argument.
    '''
    if args.fasta and args.bam:
        print 'helo'

    # Check if args.fasta is present and then call fetch_fasta with
    # args.fasta as argument.
    elif args.fasta:
        fetch_fasta(args.fasta)


    # Check if argument 'args.bam' is present and open it op for reading.
    elif args.bam:
        fetch_bam(args.bam)

# Makes sure main() is only run when this script is called from
# from itself.
if __name__ == "__main__":
    main()
