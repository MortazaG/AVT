#!/usr/bin/env python

# Licence:

# Usage:

import sys
import pysam
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

    # Ask user if they wan't to fetch information from FASTA file
    opt1 = raw_input('Fetch information from %s? [y/n]: ' % args.fasta).lower()

    if opt1 == 'y':

        with open('results/' + args.fasta + '.txt', 'w') as save_fasta:

            # Use SeqIO from BipPython to parse fasta file.
            for fasta_record in SeqIO.parse(ff, 'fasta'):
                gc_content = SeqUtils.GC(fasta_record.seq)

                # Write information from FASTA file to txt file.
                save_fasta.write('ID: %s\r\n' % fasta_record.id)
                save_fasta.write('GC content: %.2f%%\r\n' % gc_content)

            print '\n-> Results from FASTA file have been written to: results/%s.txt\n'\
                    % args.fasta

def fetch_bam(bf):

    '''
    Fetch and print necessary information from a BAM file.
    BAM filehandle is received through the bf argument.
    '''



    # Ask user if they wan't to fetch information from BAM file
    opt1 = raw_input('Fetch information from %s? [y/n]: ' % args.bam).lower()

    samfile = pysam.AlignmentFile(args.bam, 'rb')

    if opt1 == 'y':

        # Check if BAM file has an index file. If False, give option
        # to create index file.
        if samfile.has_index() == False:
            print '\nNo BAM index file could be found for %s' % (args.bam)

            bam_bai = raw_input('Produce index file now? [y/n]: ').lower()

            if bam_bai == 'y':
                pysam.index(args.bam)
                print '\nIndex file was successfully produced as: %s.bai\n'\
                        % args.bam
                exit()

            elif bam_bai == 'n':
                print '\nThis program can\'t be run without a BAM index file.\n'

            else:
                print '\nSomething wen\'t wrong...\n'

        else:

            # Calculate the total length of the reference sequence.
            # Achieved by summing up the lengths of individual references.
            ref_size = [length for length in samfile.lengths]
            ref_size = sum(ref_size)

            # Fetch amount of reads that are available.
            tot_reads = samfile.count()
            map_reads = samfile.mapped
            unmap_reads = samfile.unmapped

            # Calculate the sum of the reads lengths
            sum_reads_lengths = [len(read.seq) for read in samfile.fetch()]
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

                print '\n-> Results from BAM file have been written to: results/%s.txt\n'\
                        % args.bam

    # Call on bam_graphs to produce necessary graphs for BAM file
    bam_graphs(samfile)

    samfile.close()

def bam_graphs(sf):

    '''
    Produce necessary graphs for BAM file.
    - Coverage per reference
    - Coverage per position for each individual reference
    '''

    import matplotlib.pyplot as plt

    # Produce 'coverage per reference' graph, following:
    # coverage = (read_length*read_count)/reference_length
    # We will be using the average read length for each reference

    # Create lists to hold reference names and lengths
    refs = sf.references
    refs_lengths = sf.lengths

    # Allow user to choose which of the graphs will be produced
    opt1 = raw_input('Produce \'Coverage per reference\' graph? [y/n]: ').lower()
    opt2 = raw_input('Produce \'Coverage per position\' graph? [y/n]: ').lower()

    if opt1 == 'y':

        # Initialize list which will contain coverage for each reference
        refs_coverage = []

        # Loop through references by index nr and calculate coverage
        for i in range(len(refs)):

            # Define sum of read length for reference, so that we can calculate
            # average read length. This is done using pysams fetch()
            sum_read_length = [len(read.seq) for read in sf.fetch(refs[i])]
            sum_read_length = sum(sum_read_length)

            # Create list to hold total reads for reference
            tot_reads = sf.count(refs[i])
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

        plt.savefig('results/graphs/' + args.bam + '_cov_ref.pdf', bbox_inches='tight')
        plt.close()

        print '\n->\'Coverage per reference graph\' produced in results/graphs/'

    if opt2 == 'y':

        # Produce 'coverage per position' graph for each individual reference
        for ref in refs:

            # Use pysam pileup to get coverage for each position in reference
            ref_pos = [p.pos for p in sf.pileup(ref)]
            ref_pos_coverage = [p.n for p in sf.pileup(ref)]

            # Plot graph using matplotlib
            plt.figure(1)
            plt.plot(ref_pos, ref_pos_coverage)

            plt.title('Coverage across reference ' + ref + '\n')
            plt.ylabel('Mean Coverage')
            plt.xlabel('Position (bp)')
            plt.grid()

            plt.savefig('results/graphs/' + ref + '_cov.pdf', bbox_inches='tight')
            plt.close()

        print '->\'Coverage per position\' graph produced in results/graphs/\n'

def main():

    '''
    Main function reads in FASTA and BAM files by calling fetch_fasta()
    and fetch_bam() functions with the user input as argument.
    '''

    # Future capability to be added
    if args.fasta and args.bam:
        print 'Sorry, can\'t do that yet...'

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
