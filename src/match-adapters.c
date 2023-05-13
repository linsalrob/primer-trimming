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
#include <zlib.h>
#include "kseq.h"
#include "version.h"
#include "compare-seqs.h"
#include "print-sequences.h"
#include "reverse-complement.h"
#include "colours.h"

KSEQ_INIT(gzFile, gzread);

void encode_primers(char* primerfile, kmer_bst_t* primers, int kmer, bool reverse, int verbose) {
	/*
	 * encode the primers in primerfile 
	 * We need to use a fixed kmer length (kmer), and the primers should not be 
	 * shorter than this, so if they are we put a warning out
	 *
	 * If bool is true, we will include the reverse complement of the primers
	 */
	if( access( primerfile, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, primerfile, ENDC);
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
		uint64_t enc = kmer_encoding(seq->seq.s, 0, kmer);
		add_primer(enc, seq->name.s, primers);
		if (reverse) {
			enc = reverse_complement(enc, kmer);
			char revname[strlen(seq->name.s)+3];
			strcpy(revname, seq->name.s);
			strcat(revname, " rc");
			if (verbose)
				fprintf(stderr, "%sAdded a rc primer: %s %s\n", GREEN, revname, ENDC);
			add_primer(enc, revname, primers);
		}
			
		if (verbose)
			fprintf(stderr, "%sEncoding %s with length %ld using k-mer %d%s\n", GREEN, seq->seq.s, seq->seq.l, kmer, ENDC);
	}
}

void search_seqfile_for_primers(char* seqfile, kmer_bst_t* primers, int kmer, int verbose) {
	/* 
	 * Search through the sequences in seqfile and see if they have the 
	 * primers in primers
	 */
	if( access( seqfile, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, seqfile, ENDC);
		return;
	}

	kseq_t *seq;

	gzFile fp = gzopen(seqfile, "r");
	if (fp == NULL) {
		fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, seqfile, ENDC);
		exit(3);
	}

	seq = kseq_init(fp);
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		uint64_t enc = kmer_encoding(seq->seq.s, 0, kmer);
		kmer_bst_t *ks = find_primer(enc, primers);
		if (ks) 
			printf("%s\t%s0\n", ks->id, seq->name.s);
		
			// printf("Found primer %s (val: %ld) in sequence %s at position 0\n", ks->id, ks->value, seq->name.s);
		for (int i=1; i<seq->seq.l - kmer; i++) {
			enc = next_kmer_encoding(seq->seq.s, i, kmer, enc);
			kmer_bst_t *ks = find_primer(enc, primers);
			if (ks) 
				printf("%s\t%s\t%d\n", ks->id, seq->name.s, i);
		
			// printf("Found primer %s in sequence %s at position %d\n", ks->id, seq->name.s, i);
		}
	}
}


void trim_seq_at_primers(char* seqfile, kmer_bst_t* primers, int kmer, char* outfile, int verbose) {
	/*
	 * look through the sequences and trim at the start of the primers
	 * The same as above, but we set seq->seq.s and seq->qual.s to '\0' and write them to a file
	 */

	if( access( seqfile, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, seqfile, ENDC);
		return;
	}

	kseq_t *seq;

	gzFile fp = gzopen(seqfile, "r");
	if (fp == NULL) {
		fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, seqfile, ENDC);
		exit(3);
	}

	seq = kseq_init(fp);
	int l;
	gzFile writer = gzopen(outfile, "w");
	if (writer == NULL) {
		fprintf(stderr, "%sCan't write to %s%s\n", RED, outfile, ENDC);
		exit(4);
	}

	while ((l = kseq_read(seq)) >= 0) {
		uint64_t enc = kmer_encoding(seq->seq.s, 0, kmer);
		kmer_bst_t *ks = find_primer(enc, primers);
		if (ks) {
			fprintf(stderr, "%sIn sequence%s we found primer %s right at the beginning!!%s\n", YELLOW, seq->name.s, ks->id, ENDC);
			fprintf(stderr, "%s AT the moment, we didn't write it, so you fastq files are no longer paired!%s\n", RED, ENDC);
		} else {
			for (int i=1; i<seq->seq.l - kmer; i++) {
				enc = next_kmer_encoding(seq->seq.s, i, kmer, enc);
				kmer_bst_t *ks = find_primer(enc, primers);
				if (ks) {
					printf("%s\t%s\t%d\n", ks->id, seq->name.s, i);
					char* newseq = strdup(seq->seq.s);
					newseq[i] = '\0';
					char* newqua = strdup(seq->qual.s);
					newqua[i] = '\0';
					gzwrite(writer, "@", 1);
					gzwrite(writer, seq->name.s, strlen(seq->name.s));
					gzwrite(writer, " ", 1);
					gzwrite(writer, seq->comment.s, strlen(seq->comment.s));
					gzwrite(writer, "\n", 1);
					gzwrite(writer, newseq, i);
					gzwrite(writer, "\n+\n", 3);
					gzwrite(writer, newqua, i);
					gzwrite(writer, "\n", 1);
					break;
				}

			}
		}
	}
	gzclose(writer);


}


