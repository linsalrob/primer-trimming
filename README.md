# PRINSEQ Primer Trimming

This is a standalone package for the identification and removal of primer and adapter sequences from your fastq files. It is written in ANSI C and designed to be fast, memory efficient, and portable. 

## What is it?
Most DNA sequencing approaches use _artificial sequences_ such as adapters, primers, linkers, or MID tags, small sequences that are added to the beginning or ends of the DNA samples before they are sequenced. Here, we call them _artificial sequences_ so we don't have to keep writing adapters, primers, linkers, or MID tags!

This code is designed to identify and remove those sequences before you perform any downstream processing on the sequences.

## How does it work?

There are two steps, and one of those is optional!
1. Identify the _artifical sequences_ in your fastq file. If you know what those are (e.g. you have a file from your sequence provider, you can supply that in step 2). If you don't know what they are, we can quickly identify them for you.
2. Remove those sequences from your fastq file. We trim the primer sequences off of the left (and/or right ends).

### Predicting primers

Starting with a fastq (or fasta) file of sequences, use `primer-predictions` to identify _artificial sequences_ at the 5' end of your reads. There are a couple of input paramters you can play with:

 - `-k` to adjust the base _k_-mer size that we start with. The default (8) is a good starting point as it will allow of for the occasional sequencing error in your primers and still catch them. Its probably better not to go above 12 or 15 on a large sequence file.
 - `-p` will print out the primers and their abundances in the sequences. This is a great sanity check to make sure what we are suggesting is real (especially with 3' sequences, see below).
  - `-f` will print the primer sequences in fasta format that can be used in `primer-trimming` see below.
  - `-m` the percentage of the sequences that a _k_-mer should be present in to be included in the search. This can be a number between 1 and 100. The default is 1% of the sequences.
  