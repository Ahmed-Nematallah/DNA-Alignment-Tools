
CC=gcc
CFLAGS=-O3 -Wall -Wextra -Wfloat-equal -Wshadow -Wformat=2 \
       -Wuninitialized -Wpointer-arith -Wcast-align \
	   -Wstrict-prototypes -Wwrite-strings -Waggregate-return \
	   -Wcast-qual -Wswitch-default -Wswitch-enum \
	   -Wunreachable-code -march=native -std=c99 -pedantic \
	   -Wunused -Winvalid-pch -Wlogical-op -Wno-overlength-strings
SRCDIR=src
OBJDIR=obj
BINDIR=bin
DEPS=$(SRCDIR)/*.h

SOURCES=$(shell find $(SRCDIR)/ -name *.c)
OBJECTS=$(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

all: main

main: $(SOURCES) $(RESOURCE_OBJECTS)
	$(CC) $(SOURCES) $(RESOURCE_OBJECTS) -o $(BINDIR)/main $(CFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	rm -rf $(BINDIR)/* $(OBJDIR)/*
