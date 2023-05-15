#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "print-sequences.h"
#include "colours.h"

int encode_base(int base) {
	/*
	 * Convert a base (A, G, C, T) to a number.
	 * Note, you can use seq[i]
	 *
	 * Encoding:
	 * 	A : 0 : 00
	 * 	C : 1 : 01
	 * 	G : 2 : 10 
	 * 	T : 3 : 11
	 *
	 */

	switch ( base ) {
		case 65 :
			// A
			return 0;
		case 97 :
			// a
			return 0;
		case 67 :
			// C
			return 1;
		case 99 :
			// c
			return 1;
		case 71 :
			// G
			return 2;
		case 103 :
			// g
			return 2;
		case 84 :
			// T
			return 3;
		case 116 :
			// t
			return 3;
		case 0:
			fprintf(stderr, "%s Error: End of string (null terminator) received%s\n", RED, ENDC);
			return 0;
		default:
			fprintf(stderr, "%s ERROR: We only encode {A,G,C,T}. Don't know '%c'%s\n", RED, (char) base, ENDC);
			return 0;
	}
}

int decode_base(int val) {
	switch (val) {
		case 0:
			return 65;
			break;
		case 1:
			return 67;
			break;
		case 2:
			return 71;
			break;
		case 3:
			return 84;
			break;
		default:
			fprintf(stderr, "%s Error: We tried to decode %d but we are only using two-bit encoding.%s\n", RED, val, ENDC);
			return 78;
	}
}



uint64_t kmer_encoding(char * seq, int start_position, int k) {
	/*
	 * Given a sequence, seq, start at start_position (0 indexed), and read k characters. 
	 * Convert that to a 64-bit int using 2-bit encoding
	 *
	 * The maximum k-mer length (k-start_position is 32)
	 */

	if ((k - start_position) > 32) {
		fprintf(stderr, "%s We can only encode k<=32 strings in 64 bits. Please reduce k %s\n", PINK, ENDC);
		exit(-1);
	}

	uint64_t enc = 0;


	for (int i=start_position; i < start_position+k; i++)
		enc = (enc << 2) + encode_base(seq[i]);

	return enc;
}

char* kmer_decoding(uint64_t enc, int k) {
	/* 
	 * convert an encoded string back to a base
	 */

	char* seq = malloc(sizeof(char *) * k);
	int posn = k;
	while (posn >= 0) {
		int l = (enc & 1);
		enc >>= 1;
		int h = (enc & 1);
		enc >>= 1;
		h = (h << 1) + l;
		seq[posn--] = (char) decode_base(h);
	}
	return seq;
}



uint64_t next_kmer_encoding(char* seq, int start_position, int k, uint64_t enc) {
	/*
	 * Given a sequence, a start position, k-mer size and a previous encoding
	 * calculate the encoding that starts at start_position but remove the base
	 * at seq[start_position-1]
	 */

	if ((start_position + k - 1) > strlen(seq)) {
		fprintf(stderr, "%s Can't calculate a k-mer beyond the end of the sequence. Start: %d, k-mer %d, sequence length %ld %s\n", RED, start_position, k, strlen(seq), ENDC);
		exit(2);
	}

	enc = enc - ((uint64_t) encode_base(seq[start_position - 1]) << (2 * (k-1)));
	enc = (enc << 2) + encode_base(seq[start_position + k - 1]);
	return enc;
}





