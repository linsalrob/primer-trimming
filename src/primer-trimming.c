#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "uthash.h"

struct my_struct {
    char key[255];                    /* key */
    //int count;
    FILE *fptr;
    UT_hash_handle hh;                /* makes this structure hashable */
};


KSEQ_INIT(gzFile, gzread)

void remove_newline(char *line){
    int new_line = strlen(line) -1;
    if (line[new_line] == '\n')
        line[new_line] = '\0';
}

int trim_left(char primers[100][255], char *seq){
	int p, offset, i, match, mismatch;

	for(p = 0; p < 100; p++){
		if(primers[p][0] != '\0'){
			for(offset = 0; offset < 20; offset++){
				i = 0;
				match = mismatch = 0;
				while((mismatch<2) & (i<16)){
					if(primers[p][i] == seq[i+offset]){
						match++;
					}else{
						mismatch++;
					}
					i++;
				}
				if(match >= 15)
					return offset+16;
			}
		}
	}	
		

	return 0;

}

int main(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	struct my_struct *s;
	char primers[100][255] = { "" };
	int i, l;
	int index_pos;
	char delim[] = ":";
	int pos = 4;


	if (argc < 3) {
		fprintf(stderr, "Usage: %s <primers.tsx> <in.seq>\n", argv[0]);
		return 1;
	}


	// read in primers
	char const* const primer_file_name = argv[1];
	FILE* primer_file_pointer = fopen(primer_file_name, "r");
	char line[256];

	i = 0;
	while (fgets(line, sizeof(line), primer_file_pointer)) {
		if(line[0] != '>'){
			remove_newline(line);
			strcpy(primers[i], line);
			i++;
		}
	}
	// read in fastq
	//int line_format;
	//line_format = 0;
	fp = gzopen(argv[2], "r");
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
		index_pos = trim_left(primers, seq->seq.s);
		//printf("%s", seq->seq.s);
		for(i=index_pos; i < strlen(seq->seq.s); i++)
			putchar(seq->seq.s[i]);
		putchar('\n');

		// QUALITY
		if (seq->qual.l != seq->seq.l) continue;
		printf("+\n");
		for(i=index_pos; i < strlen(seq->seq.s); i++)
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
