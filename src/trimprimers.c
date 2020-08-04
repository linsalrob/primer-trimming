#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "kseq.h"
#include "version.h"
#include "trimprimers.h"

//#include "uthash.h"
//struct my_struct {
//    char key[255];                    /* key */
//    //int count;
//    FILE *fptr;
//    UT_hash_handle hh;                /* makes this structure hashable */
//};


KSEQ_INIT(gzFile, gzread)

#define len(x) (int)strlen(x)
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(x, y) (((x) > (y)) ? (x) : (y))

void remove_newline(char *line){
	int new_line = len(line) - 1;
	
	if(line[new_line] == '\n')
		line[new_line] = '\0';
}

int trim_left(char** primers, char *seq){
	int p, offsetS, offsetP, i, mismatch;

    for(p=0; primers[p] != NULL; p++){
		for(offsetS = 0; offsetS < 20; offsetS++){
			for(offsetP = 0; offsetP <= len(primers[p])-11; offsetP++){
				i = mismatch = 0;
				while((mismatch < 2) & (i+offsetP <= len(primers[p]))){
					if(primers[p][i+offsetP] != seq[i+offsetS])
						mismatch++;
					i++;
				}
				// this is because the last bases cannot be a mismatch
				if(primers[p][i+offsetP] == seq[i+offsetS])
					i++;
				if(i >= 11)
					return (i+offsetS-1);
			}
		}
	}	
	return 0;
}

int trim_right(char** primers, char *seq, int indexL){
	int p, offsetS, offsetP, i, mismatch;

    for(p=0; primers[p] != NULL; p++){
		for(offsetS = indexL; offsetS <= len(seq)-11; offsetS++){
			for(offsetP = 0; offsetP <= len(primers[p])-11; offsetP++){
				i = mismatch = 0;
				while((mismatch < 2) & (i+offsetP <= len(primers[p]))){
					if(primers[p][i+offsetP] != seq[i+offsetS])
						mismatch++;
					i++;
				}
				// this is because the last bases cannot be a mismatch
				if(primers[p][i+offsetP] != seq[i+offsetS])
					i--;
				// this is because the first bases cannot be a mismatch
				if(primers[p][offsetP] != seq[offsetS])
					i--;
				if(i >= 11)
					return (offsetS);
			}
		}
	}	
	return len(seq);
}

int trim_poly(char *seq, int n){
	int i = len(seq);
	char c = seq[i-1];

	while(seq[--i] == c) {}

	if(i < len(seq)-5)
		return i+1;
	else
		return len(seq);
}

char** load_primers(char *filename){
	int i, n;
	char line[256];
	FILE* fp = fopen(filename, "r");

	if(fp ==0){
		printf("Error opening file: %s\n", filename);
		exit(EXIT_FAILURE);
	}
		
	n = 0;
	while (fgets(line, sizeof(line), fp)) {
		if( (len(line) > 0) & (line[0] != '>') )
			n++;
	}
	char ** primers = malloc(++n * sizeof(char*));
	rewind(fp);

	i = 0;
	while (fgets(line, sizeof(line), fp)) {
		if(line[0] != '>'){
			remove_newline(line);
			primers[i] = malloc(255 * sizeof(char));
			strcpy(primers[i], line);
			i++;
		}
	}
	primers[i] = NULL;
	fclose(fp);
	
	return primers;
}

int trim_primers(char * infile, char **primersL, char **primersR) {
	gzFile fp;
	kseq_t *seq;
	//struct my_struct *s;
	int i, l;
	int indexL, indexR1, indexR2;


	// FASTQ
	//int line_format;
	//line_format = 0;
	fp = gzopen(infile, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		indexL = 0;
		indexR1 = indexR2 = len(seq->seq.s);
		// HEADER
		printf("%c%s", seq->qual.l == seq->seq.l? '@' : '>', seq->name.s);
		if (seq->comment.l) printf(" %s", seq->comment.s);
		putchar('\n');
		/*
		for (i = 0; i < l; ++i) {
			fputc(seq->seq.s[i], s->fptr);
		}
		*/
		// SEQUENCE
		if(primersL != NULL)
			indexL = trim_left(primersL, seq->seq.s);
		if(primersR != NULL)
			indexR1 = trim_right(primersR, seq->seq.s, indexL);
		if(1)
			indexR2 = trim_poly(seq->seq.s, 5);
		//printf("%s", seq->seq.s);
		for(i=indexL; i < min(indexR1, indexR2); i++)
			putchar(seq->seq.s[i]);
		putchar('\n');

		// QUALITY
		if (seq->qual.l != seq->seq.l) continue;
		printf("+\n");
		for(i=indexL; i < min(indexR1, indexR2); i++)
			putchar(seq->qual.s[i]);
		//printf("%s", seq->qual.s);
		/*
		for (i = 0; i < l; ++i) {
			if (i == 0 || (line_format > 0 && i % line_format == 0)) fputc('\n', s->fptr);
			fputc(seq->qual.s[i], s->fptr);
		}
		*/
		putchar('\n');
	}
	kseq_destroy(seq);
	gzclose(fp);

	return 0;
}