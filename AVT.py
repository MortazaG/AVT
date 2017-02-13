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
            sum_bases = 0
            tot_bases = 0
            for c in samfile.pileup():
                sum_bases += c.n
                tot_bases += 1

            print '\nAverage coverage/base: %.4f\n' % (float(sum_bases) / tot_bases)

        samfile.close()

if __name__ == "__main__":
    main()
