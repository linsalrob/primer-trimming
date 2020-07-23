#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include "kseq.h"
#include "uthash.h"

struct my_struct {
    char key[255];                    /* key */
    //int count;
    FILE *fptr;
    UT_hash_handle hh;                /* makes this structure hashable */
};


KSEQ_INIT(gzFile, gzread)

#define len(x) (int)strlen(x)

void remove_newline(char *line){
	int new_line = len(line) - 1;
	
	if (line[new_line] == '\n')
		line[new_line] = '\0';
}

int trim_left(char** primers, char *seq){
	int p, offsetS, offsetP, i, match, mismatch;

	for(p = 0; p < 100; p++){
		if(primers[p][0] != '\0'){
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
	}	
	return 0;
}

int trim_right(char** primers, char *seq, int indexL){
	int p, offsetS, offsetP, i, match, mismatch;

	if(0){ //for(p = 0; p < 100; p++){
		if(primers[p][0] != '\0'){
			for(offsetS = indexL; offsetS < len(seq)-11; offsetS++){
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
					if(i >= 11){
						printf("%s %i %i %i\n", primers[p], offsetS, offsetP, i);
						return (offsetS);
					}
				}
			}
		}
	}	
	return len(seq);
}

int trim_poly(char *seq, int n){
	int i;
		
	//for(i = 0; i < len(seq); i++){

	//}
	return len(seq);
}

char** load_primers(char *filename){
	int i;
	char line[256];
	FILE* fp = fopen(filename, "r");
	char ** primers = malloc(100 * sizeof(char*));
	//char primers[100][255] = { "" };

	i = 0;
	while (fgets(line, sizeof(line), fp)) {
		if(line[0] != '>'){
			remove_newline(line);
			primers[i] = malloc(255 * sizeof(char));
			strcpy(primers[i], line);
			i++;
		}
	}
	fclose(fp);
	
	return primers;
}

void print_usage() {
    printf("Usage: primer-trimming --left_primers PRIMER_FILE1 --right_primers PRIMER_FILE2 INFILE\n\n");
    printf("Primer trimming explanation...\n\n");
}

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	struct my_struct *s;
	int i, l;
	int indexL, indexR;


	// COMMAND LINE OPTIONS
	char infile[255];
	char *pfilenameL = "";
	char *pfilenameR = "";
	int opt = 0;
	static struct option long_options[] = {
	  {"left_primers",  required_argument, 0, 'l'},
	  {"right_primers",  required_argument, 0, 'r'},
	  {0, 0, 0, 0}
	};
	int option_index = 0;
	while ((opt = getopt_long(argc, argv, "l:r:", long_options, &option_index )) != -1) {
		switch (opt) {
			case 'l' : pfilenameL = optarg; 
				break;
			case 'r' : pfilenameR = optarg;
				break;
			default: print_usage(); 
				exit(EXIT_FAILURE);
		}
	}
	/* remaining command line arguments (not options). */
	if (optind < argc){
		while (optind < argc){
			strcpy(infile, argv[optind++]);
			break;
		}
 	}
	if (pfilenameL[0] == '\0' | pfilenameR[0] == '\0'){
		print_usage();
		exit(EXIT_FAILURE);
	}

	
	// PRIMERS
	//char primersL[100][255] = { "" };
	//char primersR[100][255] = { "" };
	char** primersL = load_primers(pfilenameL);
	char** primersR = load_primers(pfilenameR);


	// FASTQ
	//int line_format;
	//line_format = 0;
	fp = gzopen(infile, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {

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
		indexL = trim_left(primersL, seq->seq.s);
		indexR = trim_right(primersR, seq->seq.s, indexL);
		//printf("%s", seq->seq.s);
		for(i=indexL; i < indexR; i++)
			putchar(seq->seq.s[i]);
		putchar('\n');

		// QUALITY
		if (seq->qual.l != seq->seq.l) continue;
		printf("+\n");
		for(i=indexL; i < indexR; i++)
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

