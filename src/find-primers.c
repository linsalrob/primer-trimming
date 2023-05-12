/*
 * Find primers whereever they are in the sequence. We hash the primers and then look through the sequence to see if we have that hash
 *
 * This should be O(n) complexity where n = length of sequence
 */


#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "kseq.h"
#include "version.h"
#include "compare-seqs.h"
#include "print-sequences.h"
#include "colours.h"
#include <float.h>

KSEQ_INIT(gzFile, gzread);

void encode_primers(char* primerfile, kmer_bst_t* primers, int kmer) {
	/*
	 * encode the primers in primerfile 
	 * We need to use a fixed kmer length (kmer), and the primers should not be 
	 * shorter than this, so if they are we put a warning out
	 */

	if( access( primerfile, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "ERROR: The file %s can not be found. Please check the file path\n", primerfile);
		return;
	}

	gzFile fp;
	kseq_t *seq;

	fp = gzopen(primerfile, "r");
	seq = kseq_init(fp);
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		if (seq->seq.l < kmer) {
			fprintf(stderr, "%sWARNING: Length of %s is only %ld, so we had to adjust k-mer size down%s\n", RED, seq->name.s, seq->seq.l, ENDC);
			kmer = seq->seq.l;
		}
		printf("Encoding %s with length %ld using k-mer %d\n", seq->seq.s, seq->seq.l, kmer);
		uint64_t enc = kmer_encoding(seq->seq.s, 0, kmer);
		add_primer(enc, seq->name.s, primers);
		printf("Added primer %s with encoding %ld\n", seq->name.s, enc);
	}
}

void search_seqfile_for_primers(char* seqfile, kmer_bst_t* primers, int kmer) {
	/* 
	 * Search through the sequences in seqfile and see if they have the 
	 * primers in primers
	 */
	if( access( seqfile, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "ERROR: The file %s can not be found. Please check the file path\n", seqfile);
		return;
	}

	gzFile fp;
	kseq_t *seq;

	fp = gzopen(seqfile, "r");
	seq = kseq_init(fp);
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		uint64_t enc = kmer_encoding(seq->seq.s, 0, kmer);
		kmer_bst_t *ks = find_primer(enc, primers);
		if (ks) 
			printf("Found primer %s (val: %ld) in sequence %s at position 0\n", ks->id, ks->value, seq->name.s);
		for (int i=1; i<seq->seq.l - kmer; i++) {
			enc = next_kmer_encoding(seq->seq.s, i, kmer, enc);
			kmer_bst_t *ks = find_primer(enc, primers);
			if (ks) 
				printf("Found primer %s in sequence %s at position %d\n", ks->id, seq->name.s, i);
		}
	}
}



int main(int argc, char *argv[]) {
	char* primerfile = "primers/mgi.fa";

	kmer_bst_t *primers;
	primers = (kmer_bst_t *) malloc(sizeof(*primers));
	int kmer = 24; // our max primer length

	encode_primers(primerfile, primers, kmer);
	char* seqfile = "fastq/913873_20180417_S_R1.sample.fastq.gz";
	search_seqfile_for_primers(seqfile, primers, kmer);
}




