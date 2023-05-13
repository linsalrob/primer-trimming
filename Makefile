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

EXEC=primer-trimming primer-basecounting primer-predictions find-primers
all: $(addprefix $(BDIR), $(EXEC))


install: primer-trimming primer-predictions
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $^ $(DESTDIR)$(PREFIX)/bin


objects := $(ODIR)primer-trimming.o $(ODIR)trimprimers.o $(ODIR)primer-predictions.o $(ODIR)predictprimers.o $(ODIR)print-sequences.o \
	$(ODIR)store-primers.o $(ODIR)seqs_to_ints.o $(ODIR)find-primers.o $(ODIR)match-adapters.o

$(objects): $(ODIR)%.o: $(SDIR)%.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) $< -o $@ $(FLAGS)

$(BDIR)primer-trimming: $(SDIR)primer-trimming.o $(SDIR)trimprimers.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

primer-trimming: $(BDIR)primer-trimming

$(BDIR)primer-basecounting: $(SDIR)primer-basecounting.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

primer-basecounting: $(BDIR)primer-basecounting

$(BDIR)compare-seqs: $(SDIR)print-sequences.o $(SDIR)compare-seqs.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

compare-seqs: $(BDIR)compare-seqs

$(BDIR)primer-predictions: $(SDIR)primer-predictions.o $(SDIR)predictprimers.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

primer-predictions: $(BDIR)primer-predictions

$(BDIR)test: $(SDIR)print-sequences.o $(SDIR)store-primers.o $(SDIR)seqs_to_ints.o $(SDIR)print-sequences.o $(SDIR)test.o $(SDIR)match-adapters.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

test: $(BDIR)test

$(BDIR)find-primers: $(ODIR)print-sequences.o $(ODIR)store-primers.o $(ODIR)seqs_to_ints.o $(ODIR)find-primers.o $(ODIR)match-adapters.o
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

find-primers: $(BDIR)find-primers

.PHONY: clean

clean:
	rm -fr bin/ obj/



find-primers2: $(SDIR)print-sequences.o $(SDIR)store-primers.o $(SDIR)seqs_to_ints.o $(SDIR)match-adapters.o $(SDIR)find-primers2.o  
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)


