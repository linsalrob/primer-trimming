
/*
 * Look for adapters in R1 and R2 files at the same time.
 *
 * Very much of this code comes from [fastq-pair](https://github.com/linsalrob/fastq-pair)
 *
 * We have two files, and each file holds sequence data. Each data element is represented by only four, and exactly four
 * lines:
 *
 * @id1/1    <- the sequence identifier always begins with an @ sign and then the identifier is the string upto the first whitespace. The identifier must be unique per sequence in a file
 * ATGATC.... <- the DNA sequence
 * +        <- a line that begins with the + sign and may also have the identifier (but doesn't need to!)
 * !%$#^@     <- an ascii representation of the "quality" of the sequence. The quality is a number between 0 and 100 and this is the chr of that number +33
 * 
 * Notes:
 * 	Originally I was going to use the seek/tell approach we used in  [fastq-pair](https://github.com/linsalrob/fastq-pair) but that requires the file to be decompressed.
 * 	However, per Heng Li (http://www.htslib.org/doc/samtools-fqidx.html) this is not going to be a great solution. 
 *
 * 	> samtools fqidx should only be used on fastq files with a small number of entries. Trying to use it on a file containing millions of short sequencing reads will produce an index that is almost as big as the original file, and searches using the index will be very slow and use a lot of memory.
 *
 *	My new idea is to read the R1 file and calculate the offsets, and store those in memory, and then read R2 and do the maths. Finally, re-read the R1 file and process it (which should be quick and
 *	just I/O bound as the processing is already done).
 *
 * Approach:
 * 	Read the R1 file and store the position of each sequence and if we find an I7left primer we remember its location (Otherwise -1)
 * 	Read the R2 file and see if we find an I5right primer (Otherwise -1)
 * 	If we find neither primer, write both sequences
 *
 * 	If the offsets for the two primers match, trim both sequences.
 *
 * 	If they do not match throw an error so we get to look at those sequences to figure out what's happening
 *
 * 	If we find one, but not the other, it probably means that there is a mistake in the sequences. We should trim to the shorter of the two and mark as a potential error.
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include "match-paired-files.h"
#include "colours.h"
#include "seqs_to_ints.h"
#include "rob_dna.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread);


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
	
	if (opt->debug)
		fprintf(stderr, "%sUsing R1 kmer of %d%s\n", GREEN, r1kmer, ENDC);

	uint64_t i7l = kmer_encoding(opt->I7left, 0, r1kmer);

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
