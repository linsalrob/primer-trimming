"""
Example code to demonstrate how to incorporate predicting primers into your Python project
"""

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
    parser.add_argument('--file', '-f', help='file', required=True)
    parser.add_argument('--kmerlen', '-k', type=int, default=10,
                        help='kmer length to seed primer searching [default: %(default)d]')
    parser.add_argument('--minpercent', '-m', type=float, default=1.0,
                        help="minimum percent of reads a kmer should be in to be considered [default: %(default)d]")
    parser.add_argument('--threeprime', '-t', action="store_true", default=False,
                        help="Predict primers at the 3' end of the sequence")
    args = parser.parse_args()

    primers = PyPrinseq.primerpredict(args.file, args.kmerlen, args.minpercent, args.threeprime)

    print(f"There are {len(primers)} primers")

    for i,p in enumerate(primers):
        print(f"Python primer {i} :  {p}")
    print("Exiting")