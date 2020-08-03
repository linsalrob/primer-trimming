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
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "kseq.h"
#include "primer-predictions.h"

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

int run(char * infile, int kmerlen, double minpercent, bool print_kmer_counts, bool print_abundance,
        bool print_short_primers, bool debug) {

    // define our hash table to hold the kmers
    struct kmercount **kchash;
    kchash = malloc(sizeof(*kchash) * table_size);

    // if we are not able to allocate the memory for this, there is no point continuing!
    if (kchash == NULL) {
        fprintf(stderr, "We cannot allocate the memory for a table size of %d. Please try a smaller value\n",
                table_size);
        exit(-1);
    }

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
        printf("Reading the sequences (first time)\n");
    while ((l = kseq_read(seq)) >= 0) {
        int posn = 0;
        numseqs++;
        while (posn + kmerlen < 20) {
            char kmer[kmerlen + 1];
            substr(seq->seq.s, kmer, posn, posn + kmerlen);
            // do we have this in our array?
            unsigned h = hash(kmer) % table_size;
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

                newk->kmer = malloc(strlen(kmer) + 1);
                strcpy(newk->kmer, kmer);
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
        printf("Converting hash to array\n");

    size_t max = 0;
    size_t min = 10000000;
    struct kmercount **kcarray;
    kcarray = malloc(sizeof(*kcarray) * n);

    int p = 0;
    for (int i = 0; i < table_size; i++) {
        struct kmercount *ptr = kchash[i];
        while (ptr != NULL) {
            kcarray[p++] = ptr;
            if (sizeof(ptr) > max)
                max = sizeof(ptr);
            if (sizeof(ptr) < min)
                min = sizeof(ptr);
            ptr = ptr->next;
        }
    }
    // remove the hash now we don't need it
    free(kchash);

    qsort(kcarray, n, sizeof(*kcarray), comparator);


    if (print_kmer_counts) {
        for (int i = 0; i < n; i++) {
            printf("Kmer: %s Count: %d\n", kcarray[i]->kmer, kcarray[i]->count);
        }
    }
    fprintf(stderr, "There are %d kmers and the most appears %d times\n", n, kcarray[0]->count);

    // now we can combine adjacent kmers into longer strings
    // start with the most abundant kmer that hasn't been used

    // this is an array of char *'s primers that we match
    // we initially set it to save 100 strings, but will keep track and realloc that if required
    char **allprimers;
    int maxprimerposition = 1000;
    allprimers = malloc(sizeof(*allprimers) * maxprimerposition);
    int allprimerposition = 0;

    bool testanother = true;
    int iteration = 0;

    if (debug)
        printf("Compressing kmers\n");

    while (testanother) {
        iteration++;
        char *primer;
        int thiscount = 0;
        testanother = false;
        for (int i = 0; i < n; i++) {
            if (!kcarray[i]->used) {
                if (((double) kcarray[i]->count/numseqs) * 100 < minpercent) {
                    kcarray[i]->used = true;
                    continue;
                }
                kcarray[i]->used = true;
                primer = malloc(strlen(kcarray[i]->kmer)+1);
                strcpy(primer, kcarray[i]->kmer);
                thiscount = kcarray[i]->count;
                testanother = true;
                break;
            }
        }

        if (!*primer)
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
                    int x = 0;
                    // char tmp[strlen(primer) + 1];
                    char tmp[strlen(primer) + 1];
                    tmp[0] = kcarray[i]->kmer[0];
                    for (int j = 0; j <= strlen(primer); j++)
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
                    int j = 0;
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

        if (!*primer)
            break;

        if (strlen(primer) > kmerlen+2) {
            if (allprimerposition == maxprimerposition) {
                // realloc all primers
                maxprimerposition *= 2;
                if (debug)
                    printf("Reallocating memory for all kmers (new size: %d)\n", maxprimerposition);
                allprimers = (char **) realloc(allprimers, sizeof(*allprimers) * maxprimerposition);
            }
            allprimers[allprimerposition] = malloc(sizeof(primer));
            strcpy(allprimers[allprimerposition++], primer);

        }
        else if (print_short_primers)
                fprintf(stderr, "Skipped potential primer %s. It is too short (only %ldbp)\n", primer, strlen(primer));
        free(primer);
    }

    if (allprimerposition == 0) {
        printf("No primers could be found. It is probably because minpercent (%f) is too high. Try adding -m 0 to the command line\n", minpercent);
        return 0;
    }


    // sort the final primers by length
    if (debug)
        printf("Sorting primers\n");
    qsort(allprimers, allprimerposition-1, sizeof(*allprimers), sort_by_length);

    if (print_abundance) {
        // we are going to re-read the sequence file. I know this means two iterations, but the alternative is
        // to store the sequences as we read them, which we could do but is a pita.
        if (debug)
            printf("Printing abundance\n");
        int counts[allprimerposition];
        for (int i = 0; i<allprimerposition; i++)
            counts[i] = 0;
        gzFile fp;
        kseq_t *seq;
        //struct my_struct *s;
        int l;

        fp = gzopen(infile, "r");
        seq = kseq_init(fp);
        int maxoccurrence = 1; // the maximum value
        int n = 0; // the number of kmers we find
        while ((l = kseq_read(seq)) >= 0) {
            for (int i=0; i<allprimerposition; i++) {
                char *offset = strstr(seq->seq.s, allprimers[i]);
                if (offset) {
                    int pos = (offset - seq->seq.s) + 1;
                    if (pos < strlen(allprimers[i]) + 10) {
                        counts[i]++;
                        break;
                    }
                    else if (debug)
                        printf("For %s in %s pos %d but maxoffset %ld\n", allprimers[i], seq->name.s, pos, strlen(allprimers[i]) + 10);
                }
            }
        }
        int total = 0;
        printf("Primer\tAbundance\n");
        for (int i=0; i < allprimerposition; i++) {
            printf("%s\t%d\n", allprimers[i], counts[i]);
            total += counts[i];
        }
        printf("\nTotal primer occurrences: %d in %d sequences (%f%%)\n", total, numseqs, (float)total/numseqs*100);

    } else {
        printf("Primers found\n");
        for (int i=0; i < allprimerposition; i++)
            printf("%s\n", allprimers[i]);

    }



    return 0;
}


unsigned hash (char *s) {
    unsigned hashval;

    for (hashval=0; *s != '\0'; s++)
        hashval = *s + 31 * hashval;
    return hashval;
}

void print_usage() {
    printf("Usage: primer-predictions [OPTIONS] Sequence File (fasta or fastq)\n");
    printf("\t-k kmer length (default 8)\n");
    printf("\t-m minimum percent of the sequences that a kmer must appear in (default: 10)\n");
    printf("\t-p print abundance of each kmer\n");
    printf("\t-c print kmer counts\n\n");
    printf("\t-s print short primer sequences that are ignored\n");
    printf("Predict the primer sequences in a fasta/fastq file\n\n");
}

int main(int argc, char *argv[]) {

    // COMMAND LINE OPTIONS
    char infile[255];
    bool print_abundance = false, print_kmer_counts = false, print_short = false, debug=false;
    int kmerlen = 8;
    double minpercent = 10;
    int opt = 0;
    static struct option long_options[] = {
            {"kmer_length",  required_argument, 0, 'k'},
            {"min_percent", required_argument, 0, 'm'},
            {"print_abundance",  no_argument, 0, 'p'},
            {"print_kmer_counts",  no_argument, 0, 'c'},
            {"print_short_primers",  no_argument, 0, 's'},
            {"debug", no_argument, 0, 'd'},
            {0, 0, 0, 0}
    };
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "k:m:pcsd", long_options, &option_index )) != -1) {
        switch (opt) {
            case 'k' :
                kmerlen = atoi(optarg);
                break;
            case 'm':
                minpercent = atof(optarg);
                break;
            case 'c' : print_kmer_counts = true;
                break;
            case 'p' : print_abundance = true;
                break;
            case 's' : print_short = true;
                break;
            case 'd': debug = true;
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

    if (!*infile) {
        print_usage();
        return 0;
    }

    fprintf(stderr, "Counting kmers in %s\n", infile);
    fprintf(stderr, "kmer length: %d\nPrint kmers: %d\nPrint abundance: %i\n", kmerlen, print_kmer_counts, print_abundance);
    fprintf(stderr, "Minimum abundance percent: %f\n", minpercent);
    fprintf(stderr, "Print short primers: %d\n\n", print_short);
    return run(infile, kmerlen, minpercent, print_kmer_counts, print_abundance, print_short, debug);
}
