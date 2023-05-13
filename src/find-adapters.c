/*
 * Find primers whereever they are in the sequence. We hash the primers and then look through the sequence to see if we have that hash
 *
 * This should be O(n) complexity where n = length of sequence
 */


#include <float.h>
#include <getopt.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "colours.h"
#include "compare-seqs.h"
#include "kseq.h"
#include "match-adapters.h"
#include "print-sequences.h"
#include "version.h"


void print_usage() {
	printf("Usage: find-primers [options] -p primers.fa -f sequences.fq.gz\n");
	printf("\tprimers.fa: fasta (or fastq) file of the primer sequences to look for\n");
	printf("\tsequences.fq.gz: fastq (probably) file of the sequences to look in\n");
	printf("\tOPTIONS:\n");
	printf("\t\t-n only include forward primer sequences, not reverse complement\n");
	printf("\t\t-v print the version and exit\n");
}

int main(int argc, char *argv[]) {
	char* primerfile = NULL;
	char* fqfile = NULL;
	int opt = 0;
        static struct option long_options[] = {
                        {"primers",    required_argument, 0, 'p'},
                        {"fastq",      required_argument, 0, 'f'},
                        {"no-reverse", no_argument,       0, 'n'},
                        {"version",    no_argument,       0, 'v'},
                        {0,            0,                 0, 0}
        };
        int option_index = 0;
	bool no_reverse = false;
	bool tab = false;
        while ((opt = getopt_long(argc, argv, "p:f:nv", long_options, &option_index)) != -1) {
                switch (opt) {
                        case 'p' :
				primerfile = strdup(optarg);
                                break;
                        case 'f' :
				fqfile = strdup(optarg);
                                break;
			case 'n' :
				no_reverse = true;
				break;
                        case 'v':
                                printf("Version: %f\n", __version__);
                                return 0;
                        default:
                                print_usage();
                                exit(EXIT_FAILURE);
                }
        }

	if (primerfile == NULL && fqfile == NULL) {
		print_usage();
		exit(EXIT_FAILURE);
	}

	kmer_bst_t *primers;
	primers = (kmer_bst_t *) malloc(sizeof(*primers));
	
	int kmer = 24; // our max primer length

	bool do_reverse = true;
	if (no_reverse)
		do_reverse = false;

	encode_primers(primerfile, primers, kmer, do_reverse, 0);
	search_seqfile_for_primers(fqfile, primers, kmer, 0);

	free(primerfile);
	free(fqfile);
	free(primers);
}




