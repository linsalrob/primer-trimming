IDIR =./include
CC=gcc
CFLAGS=-g -Wall -O2 -Wno-return-type -Wno-unused-variable -Wno-unused-function -I$(IDIR)
LFLAGS= -lz -lm -g

ODIR=./src/
SDIR=./src/

LIBS=-lm


#PREFIX is environment variable, but if it is not set, then set default value
ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

all: primer-trimming primer-basecounting primer-predictions find-primers

install: primer-trimming primer-predictions
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $^ $(DESTDIR)$(PREFIX)/bin


objects = $(SDIR)primer-trimming.o $(SDIR)trimprimers.o $(SDIR)primer-predictions.o $(SDIR)predictprimers.o $(DIR)find-primers.o
$(objects): %.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@ $(FLAGS)

primer-trimming: $(SDIR)primer-trimming.c $(SDIR)trimprimers.c
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

primer-basecounting: $(SDIR)primer-basecounting.c
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

compare-seqs: $(SDIR)print-sequences.c $(SDIR)compare-seqs.c
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

primer-predictions: $(SDIR)primer-predictions.c $(SDIR)predictprimers.c
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

test: $(SDIR)print-sequences.c $(SDIR)store-primers.c $(SDIR)seqs_to_ints.c $(SDIR)print-sequences.c $(SDIR)test.c
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

find-primers: $(SDIR)print-sequences.c $(SDIR)store-primers.c $(SDIR)seqs_to_ints.c $(SDIR)print-sequences.c $(SDIR)find-primers.c
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

.PHONY: clean

clean:
	rm -f primer-trimming primer-basecounting primer-predictions src/*.o

