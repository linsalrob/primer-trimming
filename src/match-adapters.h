#ifndef MATCH_ADAPTERS_H
#define MATCH_ADAPTERS_H

/*
 * Read a primer file and convert it to encoded k-mers
 */
void encode_primers(char*, kmer_bst_t*, int, bool, int);

/*
 * Look through a fastq file and print sequences that match
 */

void search_seqfile_for_primers(char*, kmer_bst_t*, int, int);

/*
 * Look through a fastq file and then trim sequences at the first primer match
 */


void trim_seq_at_primers(char*, kmer_bst_t*, int, char*, int);

#endif
