"""
Example code to demonstrate how to incorporate trimming primers into your Python project
"""

import sys
import argparse
import PyPrinseq

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('--fastq', '-f', help='fastq file to trim primers from', required=True)
    parser.add_argument('--leftprimers', '-l', help='fasta file of left primers')
    parser.add_argument('--rightprimers', '-r', help='fasta file of right primers')
    args = parser.parse_args()

    if not args.leftprimers and not args.rightprimers:
        sys.stderr.write("ERROR: Either left or right primers must be specified. Try -h for more options\n")
        sys.exit(1)

    PyPrinseq.primertrimming(args.fastqg, args.leftprimers, args.rightprimers)
