# primer-trimming

This is completely experimental code (well, no code yet) to develop adapter primer trimming in C++.

I suggest C++ so that code can be (??) included in prinseq++ easily

The trimming is done in several steps, and we can either keep these as separate steps, write code that works the same for each step, or combine the steps. I suspect it will be easier to write one piece of code and then have it perform each step with the primer files provided as options.

For each step there is an input fastq file (e.g. `step_1_input.fq.gz`) with the sequences that should be trimmed, and an output file (e.g. `step_1_output.fq.gz`) with the output from the trimming. We may not wish to make our output exactly the same as the this output if we can justify why not!

Currently the trimming process:
- Allows a [Hamming distance](https://en.wikipedia.org/wiki/Hamming_distance) of 1
- Allows a shorter primer on the ends down to 11bp. (e.g. we could take the last 11 bp of a 16bp primer and ask if it is within 16bp of the start and then if it is > 11bp see if the others match)
- Allows the primer to be within the first/ last 20bp of the sequence



In the [primers](primers) directory there are fasta files for all the primer sequences.

### Step 1.

#### Remove leftmost primerB. Not the reverse complements

This removes the primers in `primerB.fa` from the sequences if the primer is within the first 20 bp. We do not check elsewhere in the sequences. Once a match is found, trim all bases to the left of this primer.
Do not check reverse complement.

Count how many times the first primer tag is in the input and how many times it is in the output
```sh
$ zcat fastq/step_1_input.fq.gz | grep -c TACCGTAGAGCTGCTA
3

$ ./primer-trimming primers/primerB.fa fastq/step_1_input.fq.gz | grep -c TACCGTAGAGCTGCTA
0
```

Check one of the reads to be sure that it is actually being output without the primer tag
```sh
$ zcat fastq/step_1_input.fq.gz | head -n24 | tail -n 4
@M02084:439:000000000-C4MJR:1:2112:23937:14886 1:N:0:ATCACGAT
ATACCGTAGAGCTGCTAAACCAGTACAAAAATCCAGTAAGACAAGTAAGAGTAGTATTAGTTGTTTTAGTACCAACTTTTATACTTTCAATATTACTAATAAATCTCTTAATGTCATCTTTTCTAAGAGTATCTCCAGTAGGATCTTCAAGTTTATAATAAGCCTTTTCAAATTCTTTAGCAGGACTCCAAGATTTATATCCATCTGGATAAGTTACTTCATAACCCATT
+
GGHHHGGGGGHHHHHHHHHHHHHHHHHHFHHHHHHHFGHHHHHHHGHHHGHFGFGHHHHHHGHHGHGHHEHHHHHGHGHHHHHGHHHHHGGHHHEHHHHHHHHHHHHHHHHHEFHFHFHHHHFHHHHHHHGHHFFHFGHHHGFHHHHHHDFHHHFEHGFHHFHHHHHHBFHFHGHHHBHGGGHFFHCC1G1GCGFHD<<DHB<D11GB1<F11FDF0DGF00DDGFA<<0

$ ./primer-trimming primers/primerB.fa fastq/step_1_input.fq.gz | head -n24 | tail -n 4
@M02084:439:000000000-C4MJR:1:2112:23937:14886 1:N:0:ATCACGAT
AACCAGTACAAAAATCCAGTAAGACAAGTAAGAGTAGTATTAGTTGTTTTAGTACCAACTTTTATACTTTCAATATTACTAATAAATCTCTTAATGTCATCTTTTCTAAGAGTATCTCCAGTAGGATCTTCAAGTTTATAATAAGCCTTTTCAAATTCTTTAGCAGGACTCCAAGATTTATATCCATCTGGATAAGTTACTTCATAACCCATT
+
HHHHHHHHHHHFHHHHHHHFGHHHHHHHGHHHGHFGFGHHHHHHGHHGHGHHEHHHHHGHGHHHHHGHHHHHGGHHHEHHHHHHHHHHHHHHHHHEFHFHFHHHHFHHHHHHHGHHFFHFGHHHGFHHHHHHDFHHHFEHGFHHFHHHHHHBFHFHGHHHBHGGGHFFHCC1G1GCGFHD<<DHB<D11GB1<F11FDF0DGF00DDGFA<<0

```

### Step 2.

#### Remove 3' read through contaminant

This removes the primers in `rc_primerB_ad6.fa`. Once a match is found trim all bases to the right of this primer/adapter pair. 
Do not check reverse complement.

### Step 3.

#### Remove primer free adapter (both orientations)

Remove the `nebnext_adapters.fa` from the sequences. Once found, remove sequences to the right.
Check both orientations, and check partial sequences (>10)

For example,

Before trimming:
```
@M02084:439:000000000-C4MJR:1:1106:3499:13771 1:N:0:ATCACGAT
GTGGAGTAAGTGGGCATAATGTTATGTCTACTACTCAACTTGGTGAAGCTAGAACAGGTAAATATCTTGTTATGAATAAATTTGCAGATAGTGAAACTAAACCTAGTAGACTTGATGTAGATTTTAATTCTCTCTATACTGATATTAAAAATAAATTAGTTGGAGAAGGATTTGGATCAGGTGTTGGAATTGGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATC
```

After trimming:
```
@M02084:439:000000000-C4MJR:1:1106:3499:13771 1:N:0:ATCACGAT
GTGGAGTAAGTGGGCATAATGTTATGTCTACTACTCAACTTGGTGAAGCTAGAACAGGTAAATATCTTGTTATGAATAAATTTGCAGATAGTGAAACTAAACCTAGTAGACTTGATGTAGATTTTAATTCTCTCTATACTGATATTAAAAATAAATTAGTTGGAGAAGGATTTGGATCAGGTGTTGGAATTGG
```

Removed string:
```
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATC
```

Index 1 in `nebnext_adapters.fa` is:
```
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
```

So the matched string is `AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATC` and the `ACGATCTCGTATGCCGTCTTCTGCTTG` stretch of the adapter is ignored.


### Step 4.

#### Remove adapter free primer (both orientations)

This does an exact match to the primers in `rc_primerB_ad6.fa`

In my (limited) experience, this appears to drop sequences less than 16bp long, and I didn't see any other matches.

### Other steps

Note that the other steps maybe beyond the scope of this project as we can do them with bowtie or something similar.

### Step 5: Vector contamination removal (PhiX + NCBI UniVecDB)

### Step 6: Host removal
