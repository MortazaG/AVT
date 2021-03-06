#!/usr/bin/env python
#-*- coding: UTF-8 -*-

# Licence:
# Copyright (C) 2017 Mortaza Ghanbarzehi, mortaza30@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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

# Create action to handle the -r/--rdist flag,
# using argparse.Action as parent class
class RdistAction(argparse.Action):

    # Override the __call__ method in order to set manual values,
    # namespace is the object containing all the flags and their resp. values
    def __call__(self,parser,namespace,values,option_string=None):

        # If no values are present, set rdist(self.dest) in namespace to
        # given values. self.dest = name of flag
        if not values:
            setattr(namespace,self.dest,[0, 1000, 50])

        # Otherwise, use the values that were given by the user
        else:
            setattr(namespace,self.dest,values)

# Initialization of argparser, with the name of the program and
# the necessary arguments.
parser = argparse.ArgumentParser(prog="avt.py", usage='%(prog)s <infile(s)> [optional argument(s)]')
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
                        help="Graph - plot GC%% against Coverage, requires both FASTA and BAM file. BAM file needs to be sorted by coordinate.")
parser.add_argument("-s", "--sorted", action="store_true",
                        help="Set this flag if you know for certain the BAM file is sorted")
parser.add_argument("--bcp", nargs='?', metavar='n', const=5, type=int,
                        help="BAM Graph - Coverage per position for the n longest sequences \
                        (Default = 5).")
parser.add_argument("-r", "--rdist", action=RdistAction, nargs='*', metavar='value', type=int,
                        help="BAM Graph - Histrogram distribution plot of read-pair distances.\
                        Range and bin size are set in order: min max bin (Default = 0 1000 50)")
parser.add_argument("--save", action="store_true",
                        help="Save any plotted graphs as a PDF file in graphs folder.")


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

    import numpy as np

    class _PointBrowser(object):

        """
        Click on a point to select and highlight it -- the data that
        generated the point will be shown in the lower axes.  Use the 'n'
        and 'p' keys to browse through the next and previous points
        """

        def __init__(self):
            self.lastind = 0

            self.selected, = ax.plot([refs_cov[0]], [refs_gc[0]], 'o', ms=6, alpha=0.8,
                                     color='red', visible=False)

        def onpick(self, event):

            if event.artist != plot:
                return True

            N = len(event.ind)
            if not N:
                return True

            dataind = event.ind[0]

            print 'Reference: %s\nCoverage: %.3f\nGC: %.3f%%\nLength: %d\n' % \
                        (np.take(refs_new, dataind), np.take(refs_cov, dataind), np.take(refs_gc, dataind), np.take(refs_lengths_new, dataind))

            self.lastind = dataind
            self.update()

        def update(self):
            if self.lastind is None:
                return

            dataind = self.lastind

            ax2.cla()
            ax2.plot()

            ax2.set_title('Summary for selected data point')
            ax2.text(0.05, 0.9, 'Reference: %s\nCoverage: %.3f\nGC: %.3f%%\nLength: %d\n' % \
                    (np.take(refs_new, dataind), np.take(refs_cov, dataind), np.take(refs_gc, dataind), np.take(refs_lengths_new, dataind)), \
                    transform=ax2.transAxes, va='top')

            self.selected.set_visible(True)
            self.selected.set_data(refs_cov[dataind], refs_gc[dataind])

            fig.canvas.draw()

    # Create lists to hold reference names and lengths
    refs = sf.references
    refs_lengths = sf.lengths

    # Initialize empty list which will contain coverage for each reference
    # and a list to hold the names of all refs with zero mapped bases.
    refs_cov = []
    zeros = []

    # Loop through references by index nr and calculate coverage
    for i in range(len(refs)):

        # Define sum of read length for reference, so that we can calculate
        # average read length. This is done using pysams fetch()
        read_lengths = [len(read.seq) for read in sf.fetch(refs[i])]
        sum_read_lengths = sum(read_lengths)

        if sf.count(refs[i]) < 1:
            zeros.append(refs[i])
            continue

        # Create list to hold total reads for reference
        tot_reads = sf.count(refs[i])
        av_read_length = sum_read_lengths / tot_reads

        # Define the length of current reference
        ref_length = refs_lengths[i]

        # Calculate coverage for reference according to given formula
        ref_cov = (tot_reads * av_read_length) / float(ref_length)
        refs_cov.append(ref_cov)

    # Initialize empty list to contain the gc-content of each reference sequence.
    refs_gc = []

    # Loop through each sequence to fetch its GC-content.
    for entry in SeqIO.parse(ff, 'fasta'):

        if entry.id not in zeros:
            refs_gc.append(SeqUtils.GC(entry.seq))

    # Initialize lists to contain only the refs that have mapped bases > 0
    refs_new = []
    refs_lengths_new = []

    for i in range(len(refs)):
        if refs[i] not in zeros:
            refs_new.append(refs[i])
            refs_lengths_new.append(refs_lengths[i])

    if args.save:
        ax = plt.subplot()

        ax.set_title('GC against Coverage')
        ax.set_xlabel('Coverage')
        ax.set_ylabel('GC%')

        plot, = ax.plot(refs_cov, refs_gc, 'o')

        plt.savefig('graphs/' + ff + '_GC-Coverage.pdf', bbox_inches='tight')
        plt.close()

    else:
        fig, (ax, ax2) = plt.subplots(2, 1)
        plt.subplots_adjust(hspace=0.45) # Set the distance between the two plots.

        ax.set_title('GC against Coverage')
        ax.set_xlabel('Coverage')
        ax.set_ylabel('GC%')
        ax2.tick_params(
            axis='both',
            which='both',
            bottom='off',
            top='off',
            labelbottom='off',
            labelleft='off',
            right='off',
            left='off'
        )
        plot, = ax.plot(refs_cov, refs_gc, 'o', picker=3)  # 3 points tolerance

        # Connect the user pick event with the class function browser.onpick().
        # User pick event is given as argument.
        browser = _PointBrowser()
        fig.canvas.mpl_connect('pick_event', browser.onpick)

        plt.show()

    sf.close()

def bam_cov_pos(filename, sf, top):

    '''
    Produce 'coverage per position' graph for each individual reference

    Receives filename, open_bam()'s' samfile and user input as argument.
    Outputs the n longest (user input) graphs as a pdf file in graphs/ folder.
    '''

    # Create variables to hold tuples of reference names and lengths
    refs = sf.references
    refs_lengths = sf.lengths

    # Create dictionary with the reference lengths as keys
    refs_dic = dict(zip(refs_lengths, refs))

    # Sort refs_lengths in ascending order and then store
    # the top lengths as a separate list.
    refs_lengths = sorted(refs_lengths)
    top_lengths = refs_lengths[-top:]

    # List comprehension for holding the reference names(values from refs_dic)
    # of the top lengths.
    top_refs = [refs_dic[item] for item in top_lengths]

    # Loop through individual references
    for ref in top_refs:

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
        plt.grid()

        if args.save:
            plt.savefig('graphs/' + filename + '_' + ref + '_coverage.pdf', bbox_inches='tight')
            plt.close()

        else:
            plt.show()

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

def bam_read_dist(sf, minimum, maximum, step):

    '''
    Determine the distance between read pairs.
    Take as input open_bam()'s samfile and output histogram plot.
    '''

    # Initialize lists for storing the read-pair distances
    distances = []
    negatives = []
    outliers = []

    for read in sf.fetch():

        # We don't want to get doubles, so we only use read1 and only reads
        # that are mapped.
        if read.is_read1 == True and read.is_unmapped == False:

            # This condition checks if read1's position is after its mate.
            # In that case we need to compensate in order to get a positive
            # value for the read-pair distance.
            if read.reference_start > read.next_reference_start:
                read_end = read.next_reference_start + 101
                mate_start = read.reference_start

            else:
                read_end = read.reference_end
                mate_start = read.next_reference_start

            # Distance is measured between the end position of read1 and
            # start position of its mate.
            distance = mate_start - read_end

            # Exclude distances below 0
            if distance < minimum:
                negatives.append(distance)

            # Exclude distances above 1000
            elif distance > maximum:
                outliers.append(distance)

            # Append all other distance values to distances.
            else:
                distances.append(distance)

    # Import numpy so that we can print out the mean and std.
    import numpy as np

    print
    print 'Average read-pair distance: %.2f' % np.mean(distances)
    print 'Standard deviation: %.2f' % np.std(distances)
    print 'Below %s: %d' % (minimum, len(negatives))
    print 'Above %s: %d' % (maximum, len(outliers))
    print

    # Generate bins for histogram
    bins = [x for x in range(0, max(distances), step)]

    plt.hist(distances, bins, rwidth=0.95, facecolor='g', alpha=0.8)

    plt.title('Distribution plot of read-pair distances')
    plt.xlabel('Distance (bp)')
    plt.ylabel('Number of reads')
    plt.xticks(bins)

    if args.save:
        plt.savefig('graphs/' + args.infile[0] + '_read-pair_dist.pdf', bbox_inches='tight')
        plt.close()

    else:
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

    try:
        # Check if BAM file has an index file. If False, give option
        # to create index file.
        if samfile.check_index() == False:
            print '\nNo BAM index file could be found for %s' % (bf)

            bam_bai = raw_input('Produce index file now? [y/n]: ').lower()

            if bam_bai == 'y':
                pysam.index(bf)

            elif bam_bai == 'n':
                print '\n[Error] Can\'t run program without BAM index file.\n'

            else:
                print '\nSomething wen\'t wrong...\n'

            exit()

    except AttributeError:
        return samfile

    else:
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

    import os.path

    # More than two arguments returns an Error and exits the program
    if len(infile) > 2:
        print '[Error] Too many arguments.'
        exit()

    filename = {}
    # Loop through infile list by index nr and check if it's a bam file.
    # All other files are considered to be fasta files.
    for i in range(len(infile)):

        # Check if the given filename(s) exist.
        if os.path.isfile(infile[i]) == False:
            print '\n[Error]: %s - Wrong filename / file does not exist!\n' % infile[i]
            exit()

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

            if not args.sorted:
                check_bam_sorted(filename['bam'])

            try:
                check_fasta(filename['fasta'])
                bamf_gc_cov(filename['fasta'], open_bam(filename['bam']))

            except KeyError:
                print '[Error] You forgot to type in either the BAM or FASTA filename!'

        # If True, call upon bam_cov_ref(), with filename and samfile as argument.
        if args.bcp:
            bam_cov_pos(filename['bam'], open_bam(filename['bam']), args.bcp)

        # If True, call on bam_multimapped with filename as argument.
        if args.multimapped:
            bam_multimapped(open_bam(filename['bam']))

        if args.rdist:
            bam_read_dist(open_bam(filename['bam']), args.rdist[0], args.rdist[1], args.rdist[2])


# Makes sure main() is only run when this script is called from
# from itself.
if __name__ == "__main__":

    import os

    # Create necessary folders if they don't already exist
    try:
        os.makedirs('graphs')

    except OSError:
        pass

    main()
