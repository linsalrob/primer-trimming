IDIR =./include
CC=gcc
CFLAGS=-g -Wall -O2 -Wno-return-type -Wno-unused-variable -Wno-unused-function -I$(IDIR)

SDIR=./src/

LIBS=-lm

all: primer-trimming primer-basecounting

primer-trimming: $(SDIR)primer-trimming.c
	$(CC) $(CFLAGS) $^ -o $@ -lz -lm -g

primer-basecounting: $(SDIR)primer-basecounting.c
	$(CC) $(CFLAGS) $^ -o $@ -lz -lm -g

clean:
	rm -f primer-trimming primer-basecounting
