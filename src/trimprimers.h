// Function definitions for primer.c so we can import them into python
// Created by redwards on 8/4/20.
//

#ifndef PRIMER_TRIMMING_PRIMER_TRIMMING_H
#define PRIMER_TRIMMING_PRIMER_TRIMMING_H

/*
 * Trim the primers
 */

int trim_primers(char * infile, char **primersL, char **primersR);

/*
 * trim left primers
 */

int trim_left(char** primers, char *seq);

/*
 * trim right primers
 */
int trim_right(char** primers, char *seq, int indexL);

/*
 * trim poly(A?) tails
 */
int trim_poly(char *seq, int n);


/*
 * Load the primer files
 */

char** load_primers(char *filename);


#endif //PRIMER_TRIMMING_PRIMER_TRIMMING_H
