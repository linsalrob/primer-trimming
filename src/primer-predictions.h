
#ifndef PRIMER_PREDICTIONS_H
#define PRIMER_PREDICTIONS_H

#include <stdbool.h>

/*
 * idloc is a struct with the kmer and the count of its occurence,
 * the id string for the sequence, and whether or not we've printed it out.
 * used is a boolean to determine if we have used this element in the final strings
 * next is a pointer to the next idloc element in the hash.
 */
struct kmercount {
    int count;
    char *kmer;
    bool used;
    struct kmercount *next;
};


/*
 * The method to do the running!
 *
 * Takes a char* for the file name of the fastq file, an int of the kmer length to use, a double for the minimum
 * percent of sequences that all reads should be in.
 * bool for fasta output for the primer sequences.
 * bool to print the kmer counts, and a bool to re-search through the sequences to list occurrences.
 * bool to print the short primer sequences, and a bool for debugging output
 */

int run(char * infile, int kmerlen, double minpercent, bool fasta_output, bool print_kmer_counts, bool print_abundance,
        bool print_short_primers, bool debug);

/*
 * calculate the hash for a fastq sequence
 *
 * This is a simple hash but widely used!
 *
 * we use an unsigned here so that the answer is > 0
 *
 * You still need to mod this on the table size
 */

unsigned hash (char *);


#endif //PRIMER_PREDICTIONS_H
