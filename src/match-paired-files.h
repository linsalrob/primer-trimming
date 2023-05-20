#ifndef SEARCH_PAIRED_FILES_H
#define SEARCH_PAIRED_FILES_H

#include <stdbool.h>

/*
 * The options that we need
 */

struct options {
	char* R1_file;
	char* R2_file;
	char* R1_output;
	char* R2_output;
	char* R1_matches;
	char* R2_matches;
	char* adjustments;
	char* I7left;
	char* I5right;
	int tablesize;
	bool debug;
};

/*
 * R1_read is a struct with the R1 name, the id string for the sequence, and whether there was a match to I7left
 * next is a pointer to the next R1_read element in the hash.
 */
struct R1_read {
	int trim;
	char *id;
	struct R1_read *next;
};


// how long should our lines be. This is a 64k buffer
#define MAXLINELEN 65536


/*
 * calculate the hash for a fastq sequence
 *
 * This is a simple hash but widely used!
 *
 * we use an unsigned here so that the answer is > 0
 *
 * You still need to mod this on the table size
 */

unsigned hash (char *s);

/* parse the R1 and R2 files and trim sequences */

void trim_pairwise(struct options *opt);

/* this version tries to get all single bp snps in the adapters */

void trim_pairwise_snps(struct options *opt);

#endif
