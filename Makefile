IDIR =./include
CC=gcc
CFLAGS=-g -Wall -O2 -Wno-return-type -Wno-unused-variable -Wno-unused-function -I$(IDIR)
LFLAGS= -lz -lm -g

ODIR=./obj/
SDIR=./src/
BDIR=./bin/

LIBS=-lm


#PREFIX is environment variable, but if it is not set, then set default value
ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

install: primer-trimming primer-predictions
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $^ $(DESTDIR)$(PREFIX)/bin


objects := $(ODIR)primer-trimming.o $(ODIR)trimprimers.o $(ODIR)primer-predictions.o $(ODIR)predictprimers.o $(ODIR)print-sequences.o \
	$(ODIR)store-primers.o $(ODIR)seqs_to_ints.o $(ODIR)find-adapters.o $(ODIR)match-adapters.o  $(ODIR)trim-adapters-anywhere.o \
	$(ODIR)rob_dna.o $(ODIR)primer-basecounting.o $(ODIR)test.o $(ODIR)filter_reads_with_n.o $(ODIR)match-paired-files.o \
	$(ODIR)match-paired-snps.o $(ODIR)search-paired-snp.o $(ODIR)search-paired-files.o

$(objects): $(ODIR)%.o: $(SDIR)%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@ $(FLAGS)

$(BDIR)primer-trimming: $(ODIR)primer-trimming.o $(ODIR)trimprimers.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

primer-trimming: $(BDIR)primer-trimming

$(BDIR)primer-basecounting: $(ODIR)primer-basecounting.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

primer-basecounting: $(BDIR)primer-basecounting

$(BDIR)compare-seqs: $(ODIR)print-sequences.o $(ODIR)compare-seqs.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

compare-seqs: $(BDIR)compare-seqs

$(BDIR)primer-predictions: $(ODIR)primer-predictions.o $(ODIR)predictprimers.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

primer-predictions: $(BDIR)primer-predictions

$(BDIR)test: $(ODIR)print-sequences.o $(ODIR)store-primers.o $(ODIR)seqs_to_ints.o $(ODIR)print-sequences.o $(ODIR)match-adapters.o $(ODIR)test.o $(ODIR)rob_dna.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

test: $(BDIR)test

$(BDIR)find-adapters: $(ODIR)print-sequences.o $(ODIR)store-primers.o $(ODIR)seqs_to_ints.o $(ODIR)match-adapters.o $(ODIR)find-adapters.o $(ODIR)rob_dna.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

find-adapters: $(BDIR)find-adapters
 
$(BDIR)trim-adapters-anywhere: $(ODIR)print-sequences.o $(ODIR)store-primers.o $(ODIR)seqs_to_ints.o $(ODIR)match-adapters.o $(ODIR)trim-adapters-anywhere.o $(ODIR)rob_dna.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

trim-adapters-anywhere: $(BDIR)trim-adapters-anywhere

$(BDIR)filter_reads_with_n: $(ODIR)rob_dna.o $(ODIR)filter_reads_with_n.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

filter_reads_with_n: $(BDIR)filter_reads_with_n


$(BDIR)search-paired-files: $(ODIR)match-paired-files.o $(ODIR)search-paired-files.o $(ODIR)seqs_to_ints.o $(ODIR)rob_dna.o 
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

search-paired-files: $(BDIR)search-paired-files

$(BDIR)search-paired-snps: $(ODIR)seqs_to_ints.o $(ODIR)rob_dna.o $(ODIR)store-primers.o $(ODIR)match-paired-snps.o $(ODIR)search-paired-snp.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

search-paired-snps: $(BDIR)search-paired-snps

EXEC=primer-trimming primer-basecounting primer-predictions find-adapters trim-adapters-anywhere  filter_reads_with_n
all: $(addprefix $(BDIR), $(EXEC))



.PHONY: clean

clean:
	rm -fr bin/ obj/



find-primers2: $(SDIR)print-sequences.o $(SDIR)store-primers.o $(SDIR)seqs_to_ints.o $(SDIR)match-adapters.o $(SDIR)find-primers2.o  
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)


