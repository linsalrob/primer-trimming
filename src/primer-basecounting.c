/*
 * A niaive primer prediction algorithm that just counts the bases in the first and last 20 positions of a sequence
 * and prints the most abundant base provided it is more than 50% of the composition.
 *
 * This is fast, and simple, and is a quick way to assess whether there are primers in a sequence that should
 * be removed. We do not try any heuristics to identify multiple primers (e.g. tags or indexes).
 */


#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include "kseq.h"




KSEQ_INIT(gzFile, gzread)

// how long a primer do we want to look for
#define kmer 20
// minimum fraction of sequences to be reported
#define cutoff 0.5

// we have five slots: 0: A; 1: G; 2: C; 3 T; 4: everything else
int leftCounts[kmer][5];
int numSeqs = 0;
int rightCounts[kmer][5];

void initialize_arrays() {
	for (int i=0; i<kmer; i++) {
		for (int j=0; j<5; j++) {
			leftCounts[i][j] = 0;
			rightCounts[i][j] = 0;
		}
	}
}


void countseq(char* seq) {
	// read left side primers
	for (int i=0; i<kmer; i++) {
		if (seq[i] == 65 || seq[i] == 97)
			leftCounts[i][0]++;
		else if (seq[i] == 71 || seq[i] == 103)
			leftCounts[i][1]++;
		else if (seq[i] == 67 || seq[i] == 99)
			leftCounts[i][2]++;
		else if (seq[i] == 84 || seq[i] == 116)
			leftCounts[i][3]++;
		else {
            leftCounts[i][4]++;
            printf("Adding %c to N's\n", seq[i]);
        }
	}

	// read right side primers
	int posn = kmer;
	for (int i=strlen(seq)-1; i > strlen(seq)-kmer; i--) {
		posn--;
		if (seq[i] == 65 || seq[i] == 97)
			rightCounts[posn][0]++;
		else if (seq[i] == 71 || seq[i] == 103)
			rightCounts[posn][1]++;
		else if (seq[i] == 67 || seq[i] == 99)
			rightCounts[posn][2]++;
		else if (seq[i] == 84 || seq[i] == 116)
			rightCounts[posn][3]++;
        else {
            rightCounts[posn][4]++;
            printf("Adding %c to N's\n", seq[i]);
        }
	}
}

void print_array(int counts[kmer][5]) {
    for (int i=0; i<kmer; i++) {
        printf("%d: ", i);
        for (int j=0; j<5; j++)
            printf(" %d", counts[i][j]);
        printf("\n");
    }
}

void print_possible_primer(int counts[kmer][5]) {

	for (int i = 0; i < kmer; i++) {
		// find the maximum element in the array
		int idx = 0;
		int maximum = counts[i][0];

		for (int j = 1; j < 5; j++) {
			if (counts[i][j] > maximum) {
				idx = j;
				maximum = counts[i][j];
			}
		}

		if ((float) maximum / numSeqs > cutoff) {
			if (idx == 0)
				putchar('A');
			else if (idx == 1)
				putchar('G');
			else if (idx == 2)
				putchar('C');
			else if (idx == 3)
				putchar('T');
			else {
                putchar('N');
            }
		} else {
			putchar('-');
		}
	}
	putchar('\n');
}


void print_usage() {
	printf("Usage: primer-counting INFILE\n\n");
	printf("Count the abundance of each base in the first and last 20 positions in a sequence file and print if their abundance is >0.25\n\n");
}

int main(int argc, char *argv[]){
	gzFile fp;
	kseq_t *seq;

	if (argc < 2) {
		print_usage();
		return 1;
	}

	initialize_arrays(); // set all the counts to 0

	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		countseq(seq->seq.s);
		numSeqs++;
	}
	kseq_destroy(seq);
	gzclose(fp);

	printf("Left primer: ");
	print_possible_primer(leftCounts);
	printf("Right primer: ");
	print_possible_primer(rightCounts);

	/*
	printf("Left counts\n");
	print_array(leftCounts);
	printf("\nRight counts\n");
	print_array(rightCounts);
    */

	return 0;
}

