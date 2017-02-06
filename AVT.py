#!/usr/bin/env python

# Licence:

# Usage:

import sys
from Bio import SeqIO

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
        for seq_record in SeqIO.parse(args.fasta, 'fasta'):
            print seq_record.id
            print seq_record.seq


if __name__ == "__main__":
    main()
