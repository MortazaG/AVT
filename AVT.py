#!/usr/bin/env python

# Licence:

# Usage:

import sys
import pysam
from Bio import SeqIO, SeqUtils

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

parser = argparse.ArgumentParser(prog="avt.py")
parser.add_argument("-f", "--fasta", help="Enter name of FASTA file")
parser.add_argument("-b", "--bam", help="Enter name of BAM file")
args = parser.parse_args()

def main():

    if args.fasta:
        for fasta_record in SeqIO.parse(args.fasta, 'fasta'):
            GC_content = SeqUtils.GC(fasta_record.seq)
            print 'ID: %s' % fasta_record.id
            print 'GC content: %.2f%%\n' % GC_content

    if args.bam:
        samfile = pysam.AlignmentFile(args.bam, 'rb')

        if samfile.has_index() == False:
            print '\nNo BAM index file could be found for %s.\
                    \nUse \'samtools index %s\' to create an index file. \
                    \nThe resulting file should be named as such: %s.bai\n'\
                    % (args.bam, args.bam, args.bam)
        else:

            ref_size = 0
            for length in samfile.lengths:
                ref_size += length

            for read in samfile.fetch():
                read_length = len(read.seq)
                break

            sum_bases = 0
            tot_bases = 0
            for c in samfile.pileup():
                sum_bases += c.n      # c.n gives us the coverage of base at pos x
                tot_bases += 1

            tot_reads = samfile.count()
            map_reads = samfile.mapped
            unmap_reads = samfile.unmapped

            av_coverage = (read_length * map_reads) / float(ref_size)

            print 'Reference size: %d' % ref_size
            print 'Total number of reads: %d' % tot_reads
            print 'Total number of mapped reads: %d' % map_reads
            print 'Total number of unmapped reads: %d' % unmap_reads
            print 'Read length: %d\n' % read_length

            print 'Average coverage: %.4f\n' % av_coverage

            print 'Average coverage/base: %.4f\n' % (float(sum_bases) / tot_bases)

        samfile.close()

if __name__ == "__main__":
    main()
