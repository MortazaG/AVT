#!/usr/bin/env python
#-*- coding: UTF-8 -*-

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
# -gcc (bamf_gc_cov) requires a sorted BAM file, if file is not sorted, you will
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
parser.add_argument("infile", nargs="+",
                        help="Fasta and/or BAM file")
parser.add_argument("-f", "--fasta", action="store_true",
                        help="Use in combination with -a, -g, -l, -n flags")
parser.add_argument("-a", "--all", action="store_true",
                        help="FASTA stats in columns - id, gc%%, length, N-count")
parser.add_argument("-g", '--gc', action="store_true",
                        help="FASTA id and GC%% in columns")
parser.add_argument("-l", "--length", action="store_true",
                        help="FASTA id and length in columns")
parser.add_argument("-n", "--countn", action="store_true",
                        help="FASTA id and N-count in columns")
parser.add_argument("-b", "--bam", action="store_true",
                        help="BAM stats - Overview")
parser.add_argument("-m", "--multimapped", action="store_true",
                        help="Calculate the percentage of multimapped reads from BAM file")
parser.add_argument("--gcc", action="store_true",
                        help="Graph - plot GC%% against Coverage. Requires both FASTA and BAM file.")
parser.add_argument("--bcp", action="store_true",
                        help="BAM Graph - Coverage per position")
parser.add_argument("-r", "--rdistance", action="store_true",
                        help="BAM Graph - Histrogram distribution plot of read-pair distances")

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
    Print necessary information from a BAM file.

    Receives filename and open_bam()'s samfile as argument.
    Outputs information to stdout.
    '''

    # Calculate the total length of the reference sequence.
    # Achieved by summing up the lengths of individual references.
    ref_size = [length for length in sf.lengths]
    ref_size = sum(ref_size)

    # Fetch amount of reads that are available.
    map_reads = sf.mapped
    unmap_reads = sf.unmapped
    tot_reads = map_reads + unmap_reads

    # Store mapped reads, singletons.
    map_singles = [read for read in sf.fetch() if read.mate_is_unmapped == True]

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
    print 'Number of reads: %d\r\n' % tot_reads,
    print 'Mapped reads: %d (%.2f%%)\r\n' % (map_reads,\
                                        (map_reads/float(tot_reads))*100),
    print 'Unmapped reads: %d (%.2f%%)\r\n' % (unmap_reads,\
                                        (unmap_reads/float(tot_reads))*100),
    print 'Mapped reads, singletons %d\r\n' % (len(map_singles)),
    print 'Average read length: %.2f\r\n' % av_read_length,
    print 'Average coverage: %.4f\r\n\r\n' % av_coverage,

    for ref, length in zip(refs, refs_lengths):
        print 'Reference: %s    Length (bp): %s\r\n' % (ref, length),

    sf.close()

def bamf_gc_cov(ff, sf):

    '''
    Plot GC-content against Coverage for each individual reference.

    Coverage is calculated as follows:
    coverage = (read_length*read_count)/reference_length
    We will be using the average read length for each reference

    Receives fastafile and open_bam()'s' samfile as argument.
    Outputs matplotlib graph.
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

    # Loop through each sequence to fetch its ID and GC-content.
    for entry in SeqIO.parse(ff, 'fasta'):
        id_.append(entry.id)
        refs_gc.append(SeqUtils.GC(entry.seq))

    def onpick(event):

        '''
        Take pick event as argument, when user picks a specific point on
        the graph, and print out its available data.
        '''

        # Numpy is needed to fetch the data for a specific point on the graph
        import numpy as np

        # Set index ind as the index nr of the picked point on the graph
        ind = event.ind

        # np.take(array, index) is used to fetch the available data for the
        # picked point on the graph.
        print 'Reference: %s\nCoverage: %.2f\nGC: %.2f%%\nLength: %d\n' % \
                (np.take(refs, ind), np.take(refs_cov, ind), np.take(refs_gc, ind), np.take(refs_lengths, ind))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(refs_cov, refs_gc, picker=0.5)   # Picker specifies the
                                                # sensitivity of the user
                                                # pick. A lower value
                                                # narows down the area of
                                                # the data point that the user
                                                # can click on.

    # Connect the user pick event with the function onpick(). User pickevent
    # is given as argument.
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

def bam_multimapped(sf):

    '''
    Take samfile as argument and calculate the percentage of multi-mapped reads.
    Output results to stdout.
    '''

    # Create empty variables to keep count of all the paired reads + singletons
    # and the amount of reads that are secondary.
    tot = 0
    multi_mapped = 0

    for read in sf.fetch():

        tot += 1

        # Check if reads is secondary. Returns False if read is primary.
        if read.is_secondary == True:
            multi_mapped += 1

    # Print the percentage of reads
    print multi_mapped / tot

def bam_read_dist(sf):

    '''
    Determine the distance between read pairs.
    Take as input open_bam()'s samfile and output histogram plot.
    '''

    distances = []
    for read in sf.fetch():

        if read.is_read1 == True and read.is_unmapped == False:
            if read.reference_start > read.next_reference_start:
                read_end = read.next_reference_start + 101
                mate_start = read.reference_start

            else:
                read_end = read.reference_end
                mate_start = read.next_reference_start

            distance = mate_start - read_end

            if distance < 0:
                distance = 0

            elif distance > 1000:
                distance = 0

            else:
                distances.append(distance)

    import numpy as np

    print 'Average read-pair distance: %d' % np.mean(distances)
    print 'Standard deviation: %d' % np.std(distances)

    bins = [x for x in range(0, max(distances), 50)]

    plt.hist(distances, bins, rwidth=0.95, facecolor='g', alpha=0.8)

    plt.xlabel('distance (bp)')
    plt.ylabel('Amount')
    plt.xticks(bins)
    plt.show()

def sort_bam(bf):

    '''
    BAM filename is received through the bf argument.

    Inform the user that BAM file
    needs to be sorted and asks if they wan't the program to sort it
    for them. Exit when done.
    '''

    choice = raw_input('Sort BAM file by coordinate now? [y/n]: ').lower()

    if choice == 'y':
        pysam.sort(bf, '-o', 'sorted_' + bf)
        print 'sorted_' + bf + ' was sucessfully created.'
        exit()
    else:
        print '[Error] Exiting, BAM file is not sorted accordingly'
        exit()

def check_bam_sorted(bf):

    '''
    BAM filename is received through the bf argument. Checks if file is sorted
    by looking for the @HD subtag SO and then the 'coordinate' value in header.
    If none of these are found, check if 'sorted' is in BAM filename.

    Call on sort_bam(bf) if all else fails.
    '''

    # Store all the headers for bamfile
    headers = pysam.view('-H', bf)

    try:
        # Check for the existence of @HD, SO tag. Raise exception if absent.
        if not 'SO' in headers:
            raise Exception

    except Exception:
        print '[Error] Missing tag: @HD, SO'
        sort_bam(bf)

    else:

        # Check if @HD, SO is set to coordinate and if 'sorted' can be found
        # in the bf, call on sort_bam otherwise.
        if 'coordinate' not in headers:

            if 'sorted' not in bf:
                print '[Error] BAM file is not sorted by coordinate'
                sort_bam(bf)

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

def check_fasta(ff):

    '''
    Takes as input fastafile and checks if it is valid. Outputs Error
    message if file is invalid, otherwise return ff.

    Check for:
    - Non IUPAC-Characters
    '''

    # Create list containing only IUPAC Characters
    iupac = ['g', 'a', 't', 'c', 'r', 'y', 'w', 's', 'm',
                'k', 'h', 'b', 'v', 'd', 'n', 'u', '.', '-']

    # Create list containing common non-IUPAC characters
    non_iupac = ['e', 'f', 'i', 'j', 'l', 'o', 'p', 'q', 'u', 'x', 'z',
                    '*', '#', '!', '%', '_', '&', '?', '@', '£', '$', '+', '=',
                    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

    # Loop through each individual sequence in fasta file
    for record in SeqIO.parse(ff, 'fasta'):

        seq = record.seq.lower()

        # Use find() to check if there exists non-IUPAC characters in seq.
        # If any are found, program will give an error and exit.
        for char in non_iupac:
            if seq.find(char) > -1:
                print '[Error] Non-IUPAC characters were found in %s' %\
                        record.id
                exit()

    return ff

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
            fasta_stats(check_fasta(filename['fasta']))

        # If True, call upon bam_stats(), with filename and samfile as argument.
        if args.bam:
            bam_stats(filename['bam'], open_bam(filename['bam']))

        # If True, call upon bamf_gc_cov(), with filename and samfile as argument.
        if args.gcc:
            check_bam_sorted(filename['bam'])
            bamf_gc_cov(filename['bam'], check_fasta(filename['fasta']), \
                                            open_bam(filename['bam']))

        # If True, call upon bam_cov_ref(), with filename and samfile as argument.
        if args.bcp:
            bam_cov_pos(filename['bam'], open_bam(filename['bam']))

        # If True, call on bam_multimapped with filename as argument.
        if args.multimapped:
            bam_multimapped(open_bam(filename['bam']))

        if args.rdistance:
            bam_read_dist(open_bam(filename['bam']))


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
