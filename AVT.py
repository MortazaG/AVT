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
parser.add_argument("infile", nargs="+", help="Fasta and/or BAM file")
parser.add_argument("-f", "--fasta", action="store_true", help="Use in combination with -a, -g, -l, -n flags")
parser.add_argument("-a", "--all", action="store_true", help="FASTA stats in columns - id, gc%%, length, N-count")
parser.add_argument("-g", '--gc', action="store_true", help="FASTA id and GC%% in columns")
parser.add_argument("-l", "--length", action="store_true", help="FASTA id and length in columns")
parser.add_argument("-n", "--countn", action="store_true", help="FASTA id and N-count in columns")
parser.add_argument("-b", "--bam", action="store_true", help="BAM - Overview")
parser.add_argument("--gcc", action="store_true", help="Graph - plot GC%% against Coverage.\
                                                        Requires both FASTA and BAM file.")
parser.add_argument("--bcp", action="store_true", help="BAM Graph - Coverage per position")

# Parse the above arguments, so that they can be used in the script.
args = parser.parse_args()

def count_gc(record):
    return '\t%.2f%%' % SeqUtils.GC(record)

def record_length(record):
    return '\t%s' % len(record)

def count_n(record):

    seq = record.lower()

    if seq.find('n') == -1:
        return '\t%d' % 0
    else:
        return '\t%d' % seq.find('n')

def fasta_stats(ff):

    '''
    Print ID, GC, Length and N-count from FASTA file.

    Receives FASTA filename through ff argument.
    Outputs information to stdout.
    '''

    # Use SeqIO from BipPython to parse fasta file.
    for record in SeqIO.parse(ff, 'fasta'):

        print '%s' % record.id,
        seq = record.seq

        if args.all:
            sys.stdout.write(count_gc(seq) + record_length(seq) + count_n(seq))

        if args.gc:
            sys.stdout.write(count_gc(seq))

        if args.length:
            sys.stdout.write(record_length(seq))

        if args.countn:
            sys.stdout.write(count_n(seq))

        print

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

    # Print extracted information from the BAM file to txt file
    print 'Reference size: %d\r\n' % ref_size,
    print 'Total number of reads: %d\r\n' % tot_reads,
    print 'Mapped reads: %d (%.2f%%)\r\n' % (map_reads,\
                                        (map_reads/float(tot_reads))*100),
    print 'Unmapped reads: %d (%.2f%%)\r\n' % (unmap_reads,\
                                        (unmap_reads/float(tot_reads))*100),
    print 'Average read length: %.2f\r\n' % av_read_length,
    print 'Average coverage: %.4f\r\n\r\n' % av_coverage,

    for ref, length in zip(refs, refs_lengths):
        print 'Reference: %s    Length (bp): %s\r\n' % (ref, length),

    sf.close()

def bamf_gc_cov(filename, ff, sf):

    '''
    Plot GC-content against Coverage.

    Coverage is calculated as follows:
    coverage = (read_length*read_count)/reference_length
    We will be using the average read length for each reference

    Receives filename, fastafile and open_bam()'s' samfile as argument.
    Outputs graphs as pdf file in results/graphs/ folder.
    '''

    # Create lists to hold reference names and lengths
    refs = sf.references
    refs_lengths = sf.lengths

    # Initialize empty list which will contain coverage for each reference
    refs_cov = []

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
        ref_cov = (tot_reads * av_read_length) / float(ref_length)
        refs_cov.append(ref_cov)

    # Initialize empty list to contain the id and gc-content of
    # each reference sequence.
    id_ = []
    refs_gc = []

    for entry in SeqIO.parse(ff, 'fasta'):
        id_.append(entry.id)
        refs_gc.append(SeqUtils.GC(entry.seq))

    import numpy as npy

    def onpick(event):
        ind = event.ind
        print 'Reference: %s\nCoverage: %.2f\nGC: %.2f%%\nLength: %d\n' % \
                (npy.take(refs, ind), npy.take(refs_cov, ind), npy.take(refs_gc, ind), npy.take(refs_lengths, ind))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    col = ax.scatter(refs_cov, refs_gc, picker=0.5)
    fig.canvas.mpl_connect('pick_event', onpick)

    plt.show()

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
        if args.gcc:
            bamf_gc_cov(filename['bam'], filename['fasta'], open_bam(filename['bam']))

        # If True, call upon bam_cov_ref(), with filename and samfile as argument.
        if args.bcp:
            bam_cov_pos(filename['bam'], open_bam(filename['bam']))


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
