#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "colours.h"
#include "compare-seqs.h"
#include "match-adapters.h"
#include "print-sequences.h"
#include "rob_dna.h"
#include "seqs_to_ints.h"

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

	encode_primers(primerfile, primers, kmer, true, 1);
	char* seqfile = "fastq/913873_20180417_S_R1.sample.fastq.gz";
	search_seqfile_for_primers(seqfile, primers, kmer, 1);
}


void test_reverse_complement() {
	char* seq1 = "ATGC";
	char* seq2 = "GCAT";
	uint64_t enc1 =  kmer_encoding(seq1, 0, 4);
	uint64_t enc2 =  kmer_encoding(seq2, 0, 4);
	if (reverse_complement(enc1, 4) == enc2) 
		printf("Success! We did a rc\n");
	else
		printf("Can't reverse complement %s\n", seq1);
	seq1 = "ATGCATCAGCTAGCATACGTACGTA";
	seq2 = "TACGTACGTATGCTAGCTGATGCAT";
	enc1 =  kmer_encoding(seq1, 0, 25);
	enc2 =  kmer_encoding(seq2, 0, 4);
	if (reverse_complement(enc1, 4) == enc2) 
		printf("Success! We did a rc\n");
	else
		printf("Can't reverse complement %s\n", seq1);

}

void test_has_n() {
	char* non = "ATGCATGCTACGATCGACT";
	char* withn = "AGCTAGCATNAGCTAGCAT";
	if (has_n(non)) 
		printf("WRONG: %s doesn't have an N but we think it does\n", non);
	else
		printf("Correct: %s does not have an N\n", non);
	if (has_n(withn)) 
		printf("Correct: %s DOES have an N\n", withn);
	else
		printf("WRONG: %s DOES have an N but we think it does not\n", withn);

}

void test_kmer_enc_dec() {
	char* seq = "TGAAATGCTACGATCGACT";
	uint64_t enc = kmer_encoding(seq, 0, 19);
	printf("for %s we got %s\n", seq, int_to_binary(enc));
	char* dseq = kmer_decoding(enc, 18);
	if (strcmp(seq, dseq) == 0)
		printf("%s%s and %s are the same%s\n", GREEN, seq, dseq, ENDC);
	else 
		printf("%sStarted with %s ended with %s and they are not the same!%s\n", RED, seq, dseq, ENDC);
}


int main(int argc, char *argv[]) {
	// test_primers();
	// test_comparisons();
	// match_adapters();
	// test_reverse_complement();
	// test_has_n();
	test_kmer_enc_dec();
}

