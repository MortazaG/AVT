#!/usr/bin/env python

# Licence:

# Usage:
# For BAM tasks, pysam looks for a BAM index file in the same folder.
# If index file is absent, the program will offer to produce one for you.
# The new filename will end with ''.bam.bai'.
#
# Certain BAM tasks require the file to be sorted by coordinate. If it is not,
# the program will offer to produce a sorted version of the file for you.
# The new filename will start with 'sorted_'. Alternatively, if the file is
# sorted even though the @HD, SO tag is missing or set to 'unsorted', you can
# manually add 'sorted_'to your filename and the script will recognize it
# as sorted.
#
# -gcr (gc_ref()) requires a sorted BAM file, if file is not sorted, you will
# get the option to do so. If you choose to sort it, the new filename will start
# with sorted_. This new filename must be called upon next time you run -gcr.
#
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
parser.add_argument('infile', nargs='+', help='Fasta and/or BAM filename')
parser.add_argument("-f", "--fasta", action='store_true', help="FASTA - Overview")
parser.add_argument("-b", "--bam", action='store_true', help="BAM - Overview")
parser.add_argument("--bcr", action='store_true', help="BAM Graph - Coverage per reference ")
parser.add_argument("--bcp", action='store_true', help="BAM Graph - Coverage per position")
parser.add_argument("--gcr", action='store_true', help="GC per reference graph.")

# Parse the above arguments, so that they can be used in the script.
args = parser.parse_args()

def fasta_stats(ff):

    '''
    Fetch and print necessary information from a FASTA file.
    FASTA filehandle is received through the ff argument.

    Receives FASTA filename through ff argument.
    Outputs information to txt file in results/
    '''

    # Open txt file to write results in
    with open('results/' + ff + '.txt', 'w') as save_fasta:

        # Use SeqIO from BipPython to parse fasta file.
        for fasta_record in SeqIO.parse(ff, 'fasta'):
            gc_content = SeqUtils.GC(fasta_record.seq)

            # Write information from FASTA file to txt file.
            save_fasta.write('ID: %s\r\n' % fasta_record.id)
            save_fasta.write('GC content: %.2f%%\r\n' % gc_content)

def bam_stats(filename, sf):

    '''
    Write necessary information from a BAM file to txt file.
    BAM filehandle is received through the sf argument.

    Receives filename and open_bam()'s' samfile as argument.
    Outputs information to txt file in results/ folder.
    '''

    # Calculate the total length of the reference sequence.
    # Achieved by summing up the lengths of individual references.
    ref_size = [length for length in sf.lengths]
    ref_size = sum(ref_size)

    # Fetch amount of reads that are available.
    map_reads = sf.mapped
    unmap_reads = sf.unmapped
    tot_reads = map_reads + unmap_reads

    # Calculate the sum of the reads lengths
    sum_reads_lengths = [len(read.seq) for read in sf.fetch()]
    sum_reads_lengths = sum(sum_reads_lengths)

    # Define average read length
    av_read_length = sum_reads_lengths / float(map_reads)

    # Calculate average coverage for the reference genome,
    # using the average read length.
    av_coverage = (av_read_length * map_reads) / float(ref_size)

    # Create lists containing reference names and lengths
    refs = sf.references
    refs_lengths = sf.lengths

    # Write extracted information from the BAM file to txt file
    with open('results/' + filename + '.txt', 'w') as save_bam:

        save_bam.write('Reference size: %d\r\n' % ref_size)
        save_bam.write('Total number of reads: %d\r\n' % tot_reads)
        save_bam.write('Mapped reads: %d (%.2f%%)\r\n' % (map_reads,\
                                            (map_reads/float(tot_reads))*100))
        save_bam.write('Unmapped reads: %d (%.2f%%)\r\n' % (unmap_reads,\
                                            (unmap_reads/float(tot_reads))*100))
        save_bam.write('Average read length: %.2f\r\n' % av_read_length)
        save_bam.write('Average coverage: %.4f\r\n\r\n' % av_coverage)

        for ref, length in zip(refs, refs_lengths):
            save_bam.write('Reference: %s; Length (bp): %s\r\n' % (ref, length))

    sf.close()

def bam_cov_ref(filename, sf):

    '''
    Produce 'coverage per reference' graph, following:
    coverage = (read_length*read_count)/reference_length
    We will be using the average read length for each reference

    Receives filename and open_bam()'s' samfile as argument.
    Outputs graphs as pdf fsfiles in results/graphs/ folder.
    '''

    # Create lists to hold reference names and lengths
    refs = sf.references
    refs_lengths = sf.lengths

    # Initialize list which will contain coverage for each reference
    refs_coverage = []

    # Loop through references by index nr and calculate coverage
    for i in range(len(refs)):

        # Define sum of read length for reference, so that we can calculate
        # average read length. This is done using pysams fetch()
        read_lengths = [len(read.seq) for read in sf.fetch(refs[i])]
        sum_read_lengths = sum(read_lengths)

        # Create list to hold total reads for reference
        tot_reads = sf.count(refs[i])
        av_read_length = sum_read_lengths / tot_reads

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

    plt.savefig('results/graphs/' + filename + '_cov_ref.pdf', bbox_inches='tight')
    plt.close()

    sf.close()

def bam_cov_pos(filename, sf):

    '''
    Produce 'coverage per position' graph for each individual reference

    Receives filename and open_bam()'s' samfile as argument.
    Outputs graph as a pdf file in results/graphs/ folder.
    '''

    # Create lists to hold reference names and lengths
    refs = sf.references

    # Loop through individual references
    for ref in refs:

        # Use pysam pileup to get coverage for each position in reference
        # pileup() uses 0-based indexing, thus we have to add 1 for each pos
        ref_pos = [p.pos+1 for p in sf.pileup(ref)]
        ref_pos_coverage = [p.n for p in sf.pileup(ref)]

        # Plot graph using matplotlib
        plt.figure()
        plt.plot(ref_pos, ref_pos_coverage)

        plt.title('Coverage across reference ' + ref + '\n')
        plt.ylabel('Mean Coverage')
        plt.xlabel('Position (bp)')
        plt.xticks(ref_pos, [])
        plt.grid()

        plt.savefig('results/graphs/' + filename + ref + '_cov_pos.pdf', bbox_inches='tight')
        plt.close()

    sf.close()

def gc_ref(ff, sf):

    '''
    Receive fasta and samfile as argument and plot gc per refernce graph.
    Outputs graph as a pdf file in results/graphs/ folder.
    '''

    # Create lists to hold reference names and lengths
    refs = sf.references

    # Initialize list which will contain gc for each reference
    refs_gc = []

    # Loop through references and calculate GC using SeqUtils
    for ref in SeqIO.parse(ff, 'fasta'):
        refs_gc.append(SeqUtils.GC(ref.seq))

    # Produce numerical values for x axis
    refs_range = [r+1 for r in range(len(refs))]

    # Plot graph using matplotlib.pyplot
    plt.figure()
    plt.plot(refs_range, refs_gc)

    plt.title('GC per reference\n')
    plt.ylabel('GC')
    plt.xlabel('Reference')
    plt.xticks(refs_range, refs, rotation='vertical')
    plt.grid()

    plt.savefig('results/graphs/' + ff + '_gc_ref.pdf', bbox_inches='tight')
    plt.close()

    sf.close()

def sort_bam(filename):

    '''
    Takes as argument the filename. Inform the user that BAM file
    needs to be sorted and asks if they wan't the program to sort it
    for them. Exit when done.
    '''

    choice = raw_input('Sort BAM file by coordinate now? [y/n]: ').lower()

    if choice == 'y':
        pysam.sort(filename, '-o', 'sorted_' + filename)
        exit()
    else:
        print '[Error] Exiting, BAM file is not sorted accordingly'
        exit()

def check_bam_sorted(filename, sf):

    '''
    Takes as arguments filename and samfile. Checks if file is sorted by
    looking for the @HD subtag SO, alternatively by checking if 'sorted'
    is in the filename. Call on sort_bam(filename) if none of these are found.
    '''

    # Store the header for samfile
    header = sf.header

    try:
        # Check for the existence of @HD, SO tag. Raise exception if absent.
        if not header['HD']['SO']:
            raise Exception

    except Exception:
        print '[Error] Missing tag: @HD, SO'
        sort_bam(filename)

    else:

        # Check if @HD, SO is set to coordinate and if 'sorted' can be found
        # in the filename, call on sort_bam otherwise.
        if header['HD']['SO'] != 'coordinate':

            if 'sorted' not in filename:
                print '[Error] BAM file is not sorted by coordinate'
                sort_bam(filename)

def open_bam(bf):

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
            print '\n[Error] Can\'t run program without BAM index file.\n'

        else:
            print '\nSomething wen\'t wrong...\n'

        exit()

    return samfile

def check_infile(infile):

    '''
    Check infile list for number of arguments and filetypes.
    Receives a list as argument, returns dictionary containing filetype
    as key and filename as value.
    '''

    # More than two arguments returns an Error and exits the program
    if len(infile) > 2:
        print '[Error] Too many arguments.'
        exit()

    filename = {}
    # Loop through infile list by index nr and check if it's a bam file.
    # All other files are considered to be fasta files.
    for i in range(len(infile)):

        if '.bam' in infile[i]:
            filename['bam'] = infile[i]
        else:
            filename['fasta'] = infile[i]

    return filename

def main():

    '''
    Main function checks positional argument for filename(s) and which of the
    optional arguments are called upon from the user. Output is the result from
    the respective positional argument.
    '''

    if args.infile:

        # Create dictionary holding the available filenames.
        # Dictionary key is  filetype and value is filename.
        filename = check_infile(args.infile)

        # If True, call upon fasta_stats() with fasta filename as argument.
        if args.fasta:
            fasta_stats(filename['fasta'])

        # If True, call upon bam_stats(), with filename and samfile as argument.
        if args.bam:
            bam_stats(filename['bam'], open_bam(filename['bam']))

        # If True, call upon bam_cov_ref(), with filename and samfile as argument.
        if args.bcr:
            bam_cov_ref(filename['bam'], open_bam(filename['bam']))

        # If True, call upon bam_cov_ref(), with filename and samfile as argument.
        if args.bcp:
            bam_cov_pos(filename['bam'], open_bam(filename['bam']))

        # Check if args.gcr is present and call on gc_refs to produce necessary
        # graphs.
        if args.gcr:

            # Check for index file
            bam_opened = open_bam(filename['bam'])

            # Check first if Bam file is sorted.
            check_bam_sorted(filename['bam'], bam_opened)

            # Two arguments are needed, first the fasta- and then bam-filename.
            gc_ref(filename['fasta'], bam_opened)

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
