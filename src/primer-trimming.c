/*
 * The main method for primer trimming. All of the work is done in trimprimers.c
 *
 */

#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "version.h"
#include "trimprimers.h"


void print_usage() {
    printf("Usage: primer-trimming --left_primers PRIMER_FILE1 --right_primers PRIMER_FILE2 INFILE\n\n");
    printf("Primer trimming explanation...\n\n");
}

int main(int argc, char *argv[]) {
	// COMMAND LINE OPTIONS
	char infile[255];
	char **primersL = NULL;
	char **primersR = NULL;
	int opt = 0;
	static struct option long_options[] = {
			{"left_primers",  required_argument, 0, 'l'},
			{"right_primers", required_argument, 0, 'r'},
			{"version",       no_argument,       0, 'v'},
			{0,               0,                 0, 0}
	};
	int option_index = 0;
	while ((opt = getopt_long(argc, argv, "l:r:v", long_options, &option_index)) != -1) {
		switch (opt) {
			case 'l' :
				primersL = load_primers(optarg);
				break;
			case 'r' :
				primersR = load_primers(optarg);
				break;
			case 'v':
				printf("Version: %f\n", __version__);
				return 0;
			default:
				print_usage();
				exit(EXIT_FAILURE);
		}
	}
	/* remaining command line arguments (not options). */
	if (optind < argc) {
		while (optind < argc) {
			strcpy(infile, argv[optind++]);
			break;
		}
	}
	if ((primersL == NULL) & (primersR == NULL)) {
		print_usage();
		exit(EXIT_FAILURE);
	}

	return trim_primers(infile, primersL, primersR);
}