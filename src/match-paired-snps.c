
/*
 * We are building on match-paired-files.c but this time we are also
 * searching through all possible SNPs but just in the I7left/I5right primers. We use our BST from
 * earlier to build all possibilities.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>
#include "match-paired-files.h"
#include "colours.h"
#include "seqs_to_ints.h"
#include "compare-seqs.h"
#include "rob_dna.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread);

void create_all_snps(char *adapter, int start, int kmer, kmer_bst_t *primers) {
	/*
	 * calculate a SNP at every position. There is pretty much no easy way to do this rather than an O(4**k) operation, I don't think?
	 */

	uint64_t enc = kmer_encoding(adapter, 0, kmer);
	add_primers(enc, "Original", primers);
	char* tmp = malloc(sizeof(char) * strlen(adapter));
	char bases[4] = {"A", "C", "G", "T"}; 
	for (int i=0; i<=kmer; i++) {
		strcpy(tmp, adapter);
		for (int j=0; j<4; j++) {
			tmp[i] = bases[j];
			uint64_t snp= kmer_encoding(tmp, 0, kmer);
			fprintf(stderr, "%s%s: %ld%s\n", BLUE, tmp, snp, ENDC);
			char* name = malloc(sizeof(char) * 10);

			add_primers(snp, "



void trim_pairwise(struct options *opt) {
	
	typedef struct COUNTS {
		int R1_seqs;
		int R2_seqs;
		int R1_found;
		int R2_found;
		int R1_adjusted;
		int R2_adjusted;
		int R1_trimmed;
		int R2_trimmed;
		int same;
	} COUNTS;

	COUNTS counts = {};

	
	int r1kmer = strlen(opt->I7left);
	if (r1kmer > 32) {
		fprintf(stderr, "%sI7left (%s) is longer than 32bp. Only using 32bp%s\n", GREEN, opt->I7left, ENDC);
		r1kmer = 32;
		char* cpy = strdup(opt->I7left);
		cpy[32] = '\0';
		if (opt->debug)
			fprintf(stderr, "%s%s\n%s%s\n", YELLOW, opt->I7left, cpy, ENDC);
		opt->I7left = cpy;
	}
	

	// malloc for 4**k possibilities for the bst
	kmer_bst_t *i7l_primers;
	fprintf(stderr, "%sTrying to malloc %d bytes%s\n", (sizeof(*i7l_primers * pow(4, r1kmer))), ENDC);
	i7l_primers = (kmer_bst_t *) malloc(sizeof(*i7l_primers * pow(4, r1kmer)));
	if (i7l_primers == NULL) {
		fprintf(stderr, "%sCANNOT malloc for primer bst%s\n", RED, ENDC);
		exit(3);
	}


	create_all_snps(opt->I7left, 0, r1kmer, i7l_primers);

	struct R1_read **reads;
	reads = malloc(sizeof(*reads) * opt->tablesize);

	if (reads == NULL) {
		fprintf(stderr, "%sERROR: We can not allocate memory for a table size of %d. Please try a smaller value for -t%s\n", RED, opt->tablesize, ENDC);
		exit(2);
	}

	// Step 1. Read the R1 file and find the matches to I7left primer
	if( access( opt->R1_file, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, opt->R1_file, ENDC);
		return;
	}
	if( access( opt->R2_file, R_OK ) == -1 ) {
		// file doesn't exist
		fprintf(stderr, "%sERROR: The file %s can not be found. Please check the file path%s\n", RED, opt->R2_file, ENDC);
		return;
	}

	gzFile fp1 = gzopen(opt->R1_file, "r");
	if (fp1 == NULL) {
		fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->R1_file, ENDC);
		exit(3);
	}
	kseq_t *seq = kseq_init(fp1);
	int l;

	FILE *match_out;
	if (opt->R1_matches)
		match_out = fopen(opt->R1_matches, "w");
	else
		match_out = stdout;
	
	bool warning_printed = false;
	while ((l = kseq_read(seq)) >= 0) {
		// find the location of i7l if it is in this sequence
		uint64_t enc = kmer_encoding(seq->seq.s, 0, r1kmer);
		counts.R1_seqs++;
		if (!warning_printed && has_n(seq->seq.s)) {
			fprintf(stderr, "%sWARNING: sequences have an N but we don't deal with them. They are encoded as A%s\n", BLUE, ENDC);
			warning_printed = true;
		}
		struct R1_read *R1read;
		R1read = (struct R1_read *) malloc(sizeof(*R1read));
		if (R1read == NULL) {
			fprintf(stderr, "Can't allocate memory for new ID pointer\n");
			return 0;
		}
		R1read->trim = -1;
		R1read->id = strdup(seq->name.s);
		R1read->next = NULL;


		if (enc == i7l) {
			fprintf(match_out, "R1\t%s\t0\n", seq->name.s);
			counts.R1_found++;
			R1read->trim = 0;
			unsigned hashval = hash(R1read->id) % opt->tablesize;
			R1read->next = reads[hashval];
			reads[hashval] = R1read;
			continue;
		} 

		for (int i=1; i<seq->seq.l - r1kmer + 1; i++) {
			enc = next_kmer_encoding(seq->seq.s, i, r1kmer, enc);
			if (enc == i7l) {
				fprintf(match_out, "R1\t%s\t%d\n", seq->name.s, i);
				counts.R1_found++;
				R1read->trim = i;
				break;
			}
		}
		unsigned hashval = hash(R1read->id) % opt->tablesize;
		R1read->next = reads[hashval];
		reads[hashval] = R1read;
	}

	// I am going to reset kseq so we have to initiate it again later
	kseq_destroy(seq);
	gzclose(fp1);
	fclose(match_out);

	// Step 2. Read the R2 file and find the locations of I5right
	int r2kmer = strlen(opt->I5right);
	uint64_t i5r = kmer_encoding(opt->I5right, 0, r2kmer);
	i5r = reverse_complement(i5r, r2kmer);
	
	// Open R2 for reading
	gzFile fp2 = gzopen(opt->R2_file, "r");
	if (fp2 == NULL) {
		fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->R2_file, ENDC);
		exit(3);
	}
	seq = kseq_init(fp2);
	
	// Open R2 for writing
	char* pipe_file = malloc(sizeof(char) * (strlen(opt->R2_output) + 10));
	strcpy(pipe_file, "gzip - > ");
	strcat(pipe_file, opt->R2_output);

	FILE *pipe = popen(pipe_file, "w");

	// open our log files
	FILE *adjust;
	if (opt->R2_matches)
		match_out = fopen(opt->R2_matches, "w");
	else
		match_out = stdout;
	if (opt->adjustments)
		adjust = fopen(opt->adjustments, "w");
	else
		adjust = stdout;

	fprintf(adjust, "R1/R2\tSeq ID\tFrom\tTo\n");

	while ((l = kseq_read(seq)) >= 0) {
		// find the location of i5r if it is in this sequence
		counts.R2_seqs++;
		uint64_t enc = kmer_encoding(seq->seq.s, 0, r2kmer);
		int trim = -1;
		if (enc == i5r) {
			fprintf(match_out, "R2\t%s\t0\n", seq->name.s);
			counts.R2_found++;
			trim = 0;
		} else {
			for (int i=1; i<seq->seq.l - r2kmer + 1; i++) {
				enc = next_kmer_encoding(seq->seq.s, i, r2kmer, enc);
				if (enc == i5r) {
					fprintf(match_out, "R2\t%s\t%d\n", seq->name.s, i);
					counts.R2_found++;
					trim = i;
				}
			}
		}
		// we either have a value or -1 for trim.
		// Now find the matching R1
		unsigned hashval = hash(seq->name.s) % opt->tablesize;
		struct R1_read *R1 = reads[hashval];
		bool matched = false;
		while (R1 != NULL) {
			if (strcmp(R1->id, seq->name.s) == 0) {
				// we found a match
				matched = true;
				if (trim == R1->trim) {
					// nothing to do, we can process both reads
					counts.same++;
					break;
				}
				if (trim == -1 && R1->trim > -1) {
					// fprintf(stderr, "%sWe only found a match for R1, not R2 at %s. We trimmed R2 anyway%s\n", YELLOW, seq->name.s, ENDC);
					fprintf(adjust, "R2\t%s\t%d\t%d\n", seq->name.s, trim, R1->trim);
					trim = R1->trim;
					counts.R2_adjusted++;
					break;
				}
				if (R1->trim == -1 && trim > -1) {
					// fprintf(stderr, "%sWe only found a match for R2, not R1 at %s. We trimmed R1 anyway%s\n", YELLOW, seq->name.s, ENDC);
					fprintf(adjust, "R1\t%s\t%d\t%d\n", R1->id, R1->trim, trim);
					R1->trim = trim;
					counts.R1_adjusted++;
					break;
				}

				fprintf(stderr, "%sWe want to trim starting at %d from R1 and %d from R2 in %s. We went with the shorter%s\n", BLUE, R1->trim, trim, seq->name.s, ENDC);
				if (trim < R1->trim) {
					fprintf(adjust, "R1\t%s\t%d\t%d\n", R1->id, R1->trim, trim);
					R1->trim = trim;
					counts.R1_adjusted++;
					break;
				}
				if (trim > R1->trim) {
					fprintf(adjust, "R2\t%s\t%d\t%d\n", seq->name.s, trim, R1->trim);
					trim = R1->trim;
					counts.R2_adjusted++;
					break;
				}
			}
			R1 = R1->next;
		}
		if (!matched) {
			fprintf(stderr, "%s We did not find an R1 that matches %s%s\n", PINK, seq->name.s, ENDC);
		}


		if (trim > -1) {
			seq->seq.s[trim] = '\0';
			seq->qual.s[trim] = '\0';
			counts.R2_trimmed++;
		}
		fprintf(pipe, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
	}

	pclose(pipe);
	fclose(adjust);
	fclose(match_out);
	kseq_destroy(seq);
	gzclose(fp2);

	// Step 3. Reread R1 and write the left reads, trimming at (strcmp(id, seq->name.s) == 0) -> trim
	// We only need to do this if we are going to write to the file.
	

	fp1 = gzopen(opt->R1_file, "r");
	if (fp1 == NULL) {
		fprintf(stderr, "%sERROR: Can not open %s%s\n", RED, opt->R1_file, ENDC);
		exit(3);
	}
	seq = kseq_init(fp1);

	// Open R2 for writing
	pipe_file = malloc(sizeof(char) * (strlen(opt->R1_output) + 10));
	strcpy(pipe_file, "gzip - > ");
	strcat(pipe_file, opt->R1_output);

	pipe = popen(pipe_file, "w");

	while ((l = kseq_read(seq)) >= 0) {
		unsigned hashval = hash(seq->name.s) % opt->tablesize;
		struct R1_read *R1 = reads[hashval];
		while (R1 != NULL) {
			if (strcmp(R1->id, seq->name.s) == 0) {
				if (R1->trim > -1) {
					seq->seq.s[R1->trim] = '\0';
					seq->qual.s[R1->trim] = '\0';
					counts.R1_trimmed++;
				}
			}
			R1 = R1->next;
		}
		fprintf(pipe, "@%s %s\n%s\n+\n%s\n", seq->name.s, seq->comment.s, seq->seq.s, seq->qual.s);
	}

	pclose(pipe);
	kseq_destroy(seq);
	gzclose(fp1);


	fprintf(stderr, "Total sequences: R1 %d R2 %d\n", counts.R1_seqs, counts.R2_seqs);
	fprintf(stderr, "Primer found: R1 %d R2 %d\n", counts.R1_found, counts.R2_found);
	fprintf(stderr, "Same Offset: %d (includes no adapter)\n", counts.same);
	fprintf(stderr, "Adjusted offset: R1 %d R2 %d\n", counts.R1_adjusted, counts.R2_adjusted);
	fprintf(stderr, "Sequences trimmed: R1 %d R2 %d\n", counts.R1_trimmed, counts.R2_trimmed);
}


unsigned hash (char *s) {
	unsigned hashval;

	for (hashval=0; *s != '\0'; s++)
		hashval = *s + 31 * hashval;
	return hashval;
}