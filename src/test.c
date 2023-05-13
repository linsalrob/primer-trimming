#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "colours.h"
#include "compare-seqs.h"
#include "match-adapters.h"
#include "print-sequences.h"


void test_primers() {
	// populate the bst
	
	kmer_bst_t *kmers;
	kmers = (kmer_bst_t *) malloc(sizeof(*kmers));
	uint64_t is[] = {10, 12, 11, 9, 5, 8, 14, 13};
	// uint64_t is[] = {5, 8, 9, 10, 11, 12, 13, 14};
	for (int i=0; i<= 7; i++) {
		char s[100];
		sprintf(s, "NAME: %ld", is[i]);
		add_primer(is[i], s, kmers);
		printf("Added: %ld\n", is[i]);
	}

	uint64_t is2[] = {10, 12, 11, 9, 5, 8, 14, 13, 17, 2, 99};
	for (int i=0; i<= 10; i++) {
		kmer_bst_t* ks = find_primer(is2[i], kmers);
		if ( ks == NULL )
			printf("NOT found: %ld\n", is2[i]);
		else
			printf("Found: %ld:id '%s'\n", is2[i], ks->id);
	}
}

void test_comparisons() {
	
	char * seq = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
	// char * seq = "CTCTCTCTCTCTCT";
	for (int i =1; i<32; i++) {
		uint64_t enc = kmer_encoding(seq, 0, i);
		char* bin = int_to_binary(enc);
		printf("%d gives: %li : %s\n", i, kmer_encoding(seq, 0, i), bin);
	}

	printf("\n\n\nNew hashes\n\n");
	seq = "CCCGCCCCTCCactgCCCCAAAAATTTT";
	uint64_t enc = kmer_encoding(seq, 0, 3);
	printf("i: 0 str: %s enc: %ld\t%s\n", substr(seq, 0, 3), enc, int_to_binary(enc));
	for (int i = 1; i<=20; i++) {
		enc = next_kmer_encoding(seq, i, 3, enc);
		char* ss = substr(seq, i, 3);
		printf("i: %d str: %s enc: %ld\t%s\n", i, ss, enc, int_to_binary(enc));
	}

}

void match_adapters() {
	fprintf(stderr, "%sLooking for adapters\n%s", GREEN, ENDC);
	char* primerfile = "primers/mgi.fa";

	kmer_bst_t *primers;
	primers = (kmer_bst_t *) malloc(sizeof(*primers));
	fprintf(stderr, "%sSize of primers: %ld Size of struct: %ld%s\n", BLUE, sizeof(*primers), sizeof(kmer_bst_t*), ENDC);
	int kmer = 24; // our max primer length

	encode_primers(primerfile, primers, kmer, 1);
	char* seqfile = "fastq/913873_20180417_S_R1.sample.fastq.gz";
	search_seqfile_for_primers(seqfile, primers, kmer, 1);
}


int main(int argc, char *argv[]) {
	// test_primers();
	// test_comparisons();
	match_adapters();
}

