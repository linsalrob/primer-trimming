/*
 * Identify the potential primers in a sequence.
 *
 * We start by counting all kmers (default size is 8, but you can change that here), and
 * then matching overlapping sequences provided they have 90% similar occurence (this should
 * allow for some mismatches).
 *
 * We then provide the option of searching the sequences one more time to find the occurrences of
 * each primer.
 */

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include "predictprimers.h"
#include "version.h"

void print_usage() {
    printf("Usage: primer-predictions [OPTIONS] Sequence File (fasta or fastq)\n");
    printf("\t-k kmer length (default 8)\n");
    printf("\t-m minimum percent of the sequences that a kmer must appear in (default: 1%%)\n");
    printf("\t-t predict adapter sequences on the 3' end of the reads\n");
    printf("\t-f fasta output of the primer sequences\n");
    printf("\t-p print abundance of each kmer\n");
    printf("\t-v print the version and exit\n");
    printf("\nYou probably don't want to see these outputs, but they are here if you do!\n\t-c print kmer counts\n");
    printf("\t-s print short primer sequences that are ignored\n");
    printf("Predict the primer sequences in a fasta/fastq file\n\n");
}

int main(int argc, char *argv[]) {

    // COMMAND LINE OPTIONS
    char infile[255];
    bool print_abundance = false, print_kmer_counts = false, print_short = false, debug=false, fasta_output=false;
    bool three_prime = false;
    int kmerlen = 8;
    double minpercent = 1;
    int opt = 0;
    static struct option long_options[] = {
            {"kmer_length",  required_argument, 0, 'k'},
            {"min_percent", required_argument, 0, 'm'},
            {"print_abundance",  no_argument, 0, 'p'},
            {"print_kmer_counts",  no_argument, 0, 'c'},
            {"print_short_primers",  no_argument, 0, 's'},
            {"fasta_output", no_argument, 0, 'f'},
            {"three_prime", no_argument, 0, 't'},
            {"debug", no_argument, 0, 'd'},
            {"version", no_argument, 0, 'v'},
            {0, 0, 0, 0}
    };
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "k:m:pcsdftv", long_options, &option_index )) != -1) {
        switch (opt) {
            case 'k' :
                kmerlen = atoi(optarg);
                break;
            case 'm':
                minpercent = atof(optarg);
                break;
            case 'c' : print_kmer_counts = true;
                break;
            case 'p' : print_abundance = true;
                break;
            case 's' : print_short = true;
                break;
            case 'f': fasta_output = true;
                break;
            case 't': three_prime = true;
                break;
            case 'd': debug = true;
                break;
            case 'v':
                printf("Version: %f\n", __version__);
                return 0;
            default: print_usage();
                exit(EXIT_FAILURE);
        }
    }

    /* remaining command line arguments (not options). */
    if (optind < argc){
        while (optind < argc){
            strcpy(infile, argv[optind++]);
            break;
        }
    }

    if (!*infile) {
        print_usage();
        return 0;
    }

    if (debug) {
        fprintf(stderr, "Counting kmers in %s\n", infile);
        fprintf(stderr, "kmer length: %d\n", kmerlen);
        fprintf(stderr, "Minimum abundance percent: %f\n", minpercent);
        fprintf(stderr, "fasta output: %d\n", fasta_output);
        fprintf(stderr, "identify adapters on the 3' end: %d\n", three_prime);
        fprintf(stderr, "Print kmer counts: %d\n", print_kmer_counts);
        fprintf(stderr, "Print abundance: %i\n", print_abundance);
        fprintf(stderr, "Print short primers: %d\n\n", print_short);
    }

    // for the results
    char **allprimers;
    allprimers = malloc(sizeof(*allprimers) * 1000); // initializing with max 1000 primer sequences.
    int allprimerposition = 0;


    int ro = predict_primers(infile, kmerlen, minpercent, fasta_output, three_prime, print_kmer_counts,
            print_abundance, print_short, debug, allprimers, &allprimerposition);

    if (!print_abundance && !fasta_output) {
        printf("Primers found\n");
        for (int i = 0; i < allprimerposition; i++)
            printf("Primer %d: %s\n", i, allprimers[i]);
    }

    free(allprimers);
    return ro;
}