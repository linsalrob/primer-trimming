# PRINSEQ Primer Trimming

This is a standalone package for the identification and removal of primer and adapter sequences from your fastq files. It is written in ANSI C and designed to be fast, memory efficient, and portable. 

## What is it?
Most DNA sequencing approaches use _artificial sequences_ such as adapters, primers, linkers, or MID tags, small sequences that are added to the beginning or ends of the DNA samples before they are sequenced. Here, we call them _artificial sequences_ so we don't have to keep writing adapters, primers, linkers, or MID tags!

This code is designed to identify and remove those sequences before you perform any downstream processing on the sequences.

## How does it work?

There are two steps, and one of those is optional!
1. Identify the _artifical sequences_ in your fastq file. If you know what those are (e.g. you have a file from your sequence provider, you can supply that in step 2). If you don't know what they are, we can quickly identify them for you.
2. Remove those sequences from your fastq file. We trim the primer sequences off of the left (and/or right ends).

## TL;DR

**Step 1**: Use `primer-predictions` to predict the primer sequences:

```bash
./primer-predictions -f sequences.fastq > primers.fasta
```

and if you want to check for 3' adapters:

```bash
./primer-predictions -k 10 -m 10 -f -t  sequences.fastq.gz > adapters.fasta
```

**Step 2**: Use `primer-trimming` to remove those primers:

```bash
./primer-trimming -l primers.fasta -r adapters.fasta sequences.fastq.gz > trimmed.fastq
```


### Predicting primers

Starting with a fastq (or fasta) file of sequences, use `primer-predictions` to identify _artificial sequences_ at the 5' end of your reads. There are a couple of input paramters you can play with:

 - `-k` to adjust the base _k_-mer size that we start with. The default (8) is a good starting point as it will allow of for the occasional sequencing error in your primers and still catch them. Its probably better not to go above 12 or 15 on a large sequence file.
 - `-p` will print out the primers and their abundances in the sequences. This is a great sanity check to make sure what we are suggesting is real (especially with 3' sequences, see below).
  - `-f` will print the primer sequences in fasta format that can be used in `primer-trimming` see below.
  - `-m` the percentage of the sequences that a _k_-mer should be present in to be included in the search. This can be a number between 1 and 100. The default is 1% of the sequences.
  - `-t` look for adapters on the 3' end of the sequences (see below).
  
 There are some other options that are largely for debuging the code, and you are free to explore them, but you will likely not need to use or change them.
 
#### 5' and 3' Searching

We can search either the 5' end of the sequences (e.g. for MID tags, linkers, primers, etc) or the 3' end of the sequences (e.g. for a reverse primer). If you are searching the 3' end of the sequences we strongly recommend printing the counts of the primers in the actual sequences (using the `-p` option) to the code, and then deciding whether those sequences need trimming off with `primer-trimming`. In our experience we mostly don't find a 3' sequence worthy of trimming, but we may still report it here. We also recommend increasing the value of `-m` as this will remove errant sequences that are not really primers.


### Trimming primers 

We can trim either 5' primers (using the `-l` or `--left_primers` option) or 3' primers (using the `-r` or `--right_primers` option). We also trim poly-A tails off of sequences by default.

The `-l` and `-r` options need a fasta format file, and so you can use the output from `primer-prediction` above directly in the trimming step here.

## Installation

You can install from source, just clone this git repo and use make:

```bash
git clone https://github.com/linsalrob/primer-trimming.git
cd primer-trimming
make all
```

Conda installation is coming soon (bug Rob about it!)