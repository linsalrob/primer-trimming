/*
 * Identify the potential primers in a sequence.
 *
 * We start by counting all kmers (default size is 8, but you can change that here), and
 * then matching overlapping sequences provided they have 90% similar occurence (this should
 * allow for some mismatches).
 *
 * We then provide the option of searching the sequences one more time to find the occurrences of
 * each primer.
 */

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include "kseq.h"
#include "predictprimers.h"
#include "version.h"

KSEQ_INIT(gzFile, gzread)

#define table_size 10000



void substr(char* seq, char * kmer, int start, int stop) {
    int p = 0;
    for (int i=start; i<stop; i++)
        kmer[p++] = seq[i];
    kmer[p] = 0;
}

/*
 * Method to compare two kmercount objects used in the qsort method
 */

int comparator(const void *p, const void *q) {
    /*
     * Left this here to enable debugging the comparator
     * struct kmercount *kc = *((struct kmercount **)q);
     * struct kmercount *kp = *((struct kmercount **)p);
     * printf("KMER: %s Count %d  KMER: %s Count %d\n", kc->kmer, kc->count, kp->kmer, kp->count);
     */
    return ((*(struct kmercount **)q)->count - (*(struct kmercount **)p)->count);
}

int sort_by_length(const void *p, const void *q) {
    /*
     * Left this here to enable debugging the comparator
     * const char *s1 = *(const char * const *)p;
     * const char *s2 = *(const char * const *)q;
     * printf("For %s strlength is %ld\n", s1, strlen(s2));
     */
    return (strlen(*(const char * const *)q) - strlen(*(const char * const *)p));
}


/*
 * The basic concept is that we don't know - a priori - how many kmers we will find, so we make a hash
 * of kmers and their counts. max is actually 4**kmerlength and with small kmers (<=10) and moderate
 * sized sequences we typically find all kmers, but this scales to larger k.
 *
 * Once we have all kmers we put them into an array so that we can sort by frequency (kcharray[i]->count)
 * so we can combine the most abundant kmers first. A hash can't be sorted.
 *
 * We use two steps for this - convert hash to array and then qsort, but we probably could have done an
 * insertion sort too. But since the data is unsorted, probably no speed up doing insertion sort.
 *
 * We iterate through the array and try and merge kmers - making a note of ones that we have used by setting
 * their boolean.
 *
 * We then have an optional step of looking back through the sequences to see if we can find those kmers. We
 * decided not to keep all the sequences in memory, but rather iterate through the file twice as this
 * would help with large sequence files.
 *
 */

int predict_primers(char * infile, int kmerlen, double minpercent, bool fasta_output, bool three_prime, bool print_kmer_counts, bool print_abundance,
        bool print_short_primers, bool debug, char **allprimers, int *allprimerposition) {

    if( access( infile, R_OK ) == -1 ) {
        // file doesn't exist
        fprintf(stderr, "ERROR: The file %s can not be found. Please check the file path\n", infile);
        return 1;
    }


    // define our hash table to hold the kmers
    struct kmercount **kchash;
    kchash = malloc(sizeof(*kchash) * table_size);

    // if we are not able to allocate the memory for this, there is no point continuing!
    if (kchash == NULL) {
        fprintf(stderr, "We cannot allocate the memory for a table size of %d. Please try a smaller value\n",
                table_size);
        exit(-1);
    }

    // initilize the table to be empty
    for (int i = 0; i<table_size; i++)
        kchash[i] = NULL;

    gzFile fp;
    kseq_t *seq;
    //struct my_struct *s;
    int l;

    fp = gzopen(infile, "r");
    seq = kseq_init(fp);
    int maxoccurrence = 1; // the maximum value
    int n = 0; // the number of kmers we find
    int numseqs = 0;
    if (debug)
        fprintf(stderr, "Reading the sequences (first time)\n");
    while ((l = kseq_read(seq)) >= 0) {
        numseqs++;
        double posn = 0;
        if (three_prime) {
            posn = strlen(seq->seq.s) - 20 - kmerlen;
        }
        bool moreseqs = true;
        while (moreseqs) {
            if (three_prime) {
                if (posn + kmerlen == strlen(seq->seq.s))
                    moreseqs = false;
            }
            else {
                if (posn > 20)
                    moreseqs = false;
            }
            char kmer[kmerlen + 1];
            substr(seq->seq.s, kmer, posn, (posn + kmerlen));
            // do we have this in our array?
            unsigned int h = hash(kmer) % table_size;
            struct kmercount *ptr = kchash[h];
            int matched = 0;

            while (ptr != NULL) {
                if (strcmp(ptr->kmer, kmer) == 0) {
                    matched = 1;
                    ptr->count++;
                    if (ptr->count > maxoccurrence) {
                        maxoccurrence = ptr->count;
                    }
                }
                ptr = ptr->next;
            }
            if (!matched) {
                struct kmercount *newk;
                newk = (struct kmercount *) malloc(sizeof(*newk));
                newk->kmer = strdup(kmer);
                //newk->kmer = malloc(strlen(kmer) + 1);
                //strcpy(newk->kmer, kmer);

                newk->count = 1;
                newk->used = false;
                // add first for the kmer to the hash
                newk->next = kchash[h];
                kchash[h] = newk;
                n++;
            }
            posn++;
        }
    }
    kseq_destroy(seq);
    gzclose(fp);

    // now we know how many kmers we have, we can convert our hash to an
    // array of elements

    if (debug)
        fprintf(stderr, "Converting hash to array\n");

    struct kmercount **kcarray;
    kcarray = malloc(sizeof(*kcarray) * n);

    int p = 0;
    for (int i = 0; i < table_size; i++) {
        struct kmercount *ptr = kchash[i];
        while (ptr != NULL) {
            kcarray[p++] = ptr;
            ptr = ptr->next;
        }
    }
    // remove the hash now we don't need it
    free(kchash);

    if (debug)
        fprintf(stderr, "Quick sorting\n");

    qsort(kcarray, n, sizeof(*kcarray), comparator);

    if (print_kmer_counts) {
        for (int i = 0; i < n; i++) {
            printf("Kmer: %s Count: %d\n", kcarray[i]->kmer, kcarray[i]->count);
        }
    }
    if (debug)
        fprintf(stderr, "There are %d kmers and the most appears %d times\n", n, kcarray[0]->count);

    // now we can combine adjacent kmers into longer strings
    // start with the most abundant kmer that hasn't been used

    // this is an array of char *'s primers that we match
    // we initially set it to save 100 strings, but will keep track and realloc that if required

    int maxprimerposition = 1000;
    //allprimers = malloc(sizeof(*allprimers) * maxprimerposition);
    //allprimers = (char **) realloc(allprimers, sizeof(*allprimers) * maxprimerposition);
    *allprimerposition = 0;

    bool testanother = true;
    int iteration = 0;

    if (debug)
        fprintf(stderr, "Compressing kmers\n");

    char *primer = NULL;
    while (testanother) {
        iteration++;
        int thiscount = 0;
        testanother = false;
        for (int i = 0; i < n; i++) {
            if (!kcarray[i]->used) {
                if (((double) kcarray[i]->count/numseqs) * 100 < minpercent) {
                    kcarray[i]->used = true;
                    continue;
                }
                kcarray[i]->used = true;
                primer = strdup(kcarray[i]->kmer);
                thiscount = kcarray[i]->count;
                testanother = true;
                break;
            }
        }


        if (debug && !primer)
            fprintf(stderr, "No primer sequence. Breaking\n");

        if (!primer)
            break;

        if (debug)
            fprintf(stderr, "Testing %s\n", primer);

        bool matched = true;
        while (matched) {
            matched = false;
            for (int i = 0; i < n; i++) {
                //printf("Iteration: %d I: %d primer: %s kmer: %s Used: %d\n", iteration, i, primer, kcarray[i]->kmer, kcarray[i]->used);
                if (kcarray[i]->used)
                    continue;
                if (((float) kcarray[i]->count / thiscount) < 0.9)
                    continue;
                // check to see if the strings match at the beginning
                char kbeg[kmerlen + 1];
                char pbeg[kmerlen + 1];
                substr(kcarray[i]->kmer, kbeg, 1, kmerlen);
                substr(primer, pbeg, 0, kmerlen - 1);

                if (strcmp(kbeg, pbeg) == 0) {
                    // char tmp[strlen(primer) + 1];
                    char tmp[strlen(primer) + 1];
                    tmp[0] = kcarray[i]->kmer[0];
                    for (unsigned long j = 0; j <= strlen(primer); j++)
                        tmp[j + 1] = primer[j];
                    primer = (char *) realloc(primer, strlen(tmp)+1);
                    strcpy(primer, tmp);
                    kcarray[i]->used = true;

                    matched = true;
                    continue;
                }

                // check for matches at the end
                char kend[kmerlen + 1];
                char pend[kmerlen + 1];
                substr(kcarray[i]->kmer, kend, 0, kmerlen - 1);
                int pstartpos = (strlen(primer) - kmerlen) + 1;
                substr(primer, pend, pstartpos, strlen(primer));

                if (strcmp(kend, pend) == 0) {
                    char tmp[strlen(primer) + 10];
                    unsigned long j = 0;
                    while (j < strlen(primer)) {
                        tmp[j] = primer[j];
                        j++;
                    }
                    tmp[j] = kcarray[i]->kmer[strlen(kcarray[i]->kmer) - 1];
                    tmp[j + 1] = 0;
                    primer = (char *) realloc(primer, strlen(tmp)+1);
                    strcpy(primer, tmp);
                    kcarray[i]->used = true;
                    matched = true;
                    continue;
                }
            }
        }

        if (debug && !primer)
            fprintf(stderr, "No primer sequence. Breaking\n");

        if (!primer)
            break;

        if ((int) strlen(primer) > kmerlen+2) {
            bool addthis = true;
            for (int i=0; i < *allprimerposition; i++) {
                if (strcmp(allprimers[i], primer) == 0)
                    addthis = false;
            }
            if (addthis) {
                if (*allprimerposition == maxprimerposition) {
                    // realloc all primers
                    maxprimerposition *= 2;
                    if (debug)
                        fprintf(stderr, "Reallocating memory for all kmers (new size: %d)\n", maxprimerposition);
                    allprimers = (char **) realloc(allprimers, sizeof(*allprimers) * maxprimerposition);
                }
                allprimers[(*allprimerposition)++] = strdup(primer);
            }
        }
        else if (print_short_primers)
                fprintf(stderr, "Skipped potential primer %s. It is too short (only %ldbp)\n", primer, strlen(primer));

    }
    free(primer);

    if (*allprimerposition == 0) {
        printf("No primers could be found. It is probably because minpercent (%f) is too high. Try adding -m 0 to the command line\n", minpercent);
        return 0;
    }

    // sort the final primers by length
    if (debug)
        fprintf(stderr, "Sorting primers\n");
    qsort(allprimers, (*allprimerposition)-1, sizeof(*allprimers), sort_by_length);

    if (print_abundance) {
        // we are going to re-read the sequence file. I know this means two iterations, but the alternative is
        // to store the sequences as we read them, which we could do but is a pita.
        if (debug)
            fprintf(stderr, "Printing abundance\n");
        int counts[*allprimerposition];
        for (int i = 0; i<*allprimerposition; i++)
            counts[i] = 0;
        gzFile fp;
        kseq_t *seq;
        //struct my_struct *s;
        int l;

        fp = gzopen(infile, "r");
        seq = kseq_init(fp);
        while ((l = kseq_read(seq)) >= 0) {
            for (int i=0; i<*allprimerposition; i++) {
                char *offset = strstr(seq->seq.s, allprimers[i]);
                if (offset) {
                    unsigned long pos = (offset - seq->seq.s) + 1;
                    if (three_prime) {
                        if (debug)
                            fprintf(stderr, "For %s in %s pos %ld but strlen %ld and maxoffset %ld\n", allprimers[i], seq->name.s, pos, strlen(seq->seq.s), strlen(allprimers[i]) + 10);
                        if (strlen(seq->seq.s) - pos < (unsigned long) (20 + kmerlen)) {
                            counts[i]++;
                            break;
                        }
                    } else {
                        if (pos < strlen(allprimers[i]) + 20) {
                            counts[i]++;
                            break;
                        }
                    }
                    if (debug)
                        fprintf(stderr, "For %s in %s pos %ld but maxoffset %ld\n", allprimers[i], seq->name.s, pos, strlen(allprimers[i]) + 10);
                }
            }
        }
        int total = 0;
        printf("Primer\tAbundance\n");
        for (int i=0; i < *allprimerposition; i++) {
            printf("%s\t%d\n", allprimers[i], counts[i]);
            total += counts[i];
        }
        printf("\nTotal primer occurrences: %d in %d sequences (%f%%)\n", total, numseqs, (float)total/numseqs*100);

    }
    if (fasta_output) {
        for (int i = 0; i < *allprimerposition; i++)
            printf(">primer_%d\n%s\n", i, allprimers[i]);
    }
    return 0;
}


unsigned int hash (char *s) {
    unsigned int hashval;

    for (hashval=0; *s != '\0'; s++)
        hashval = *s + 31 * hashval;
    return hashval;
}
