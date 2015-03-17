PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

BOOSTDIR=${HOME}/homebrew

CXXFLAGS=-Wall -g -std=c++11 -stdlib=libc++
CPPFLAGS=-I$(BOOSTDIR)/include
LDLIBS=-L$(BOOSTDIR)/lib -lm -lboost_program_options

PROGRAMS=pangwas
CLASSES=sample.o kmer.o
OBJECTS=pangwasMain.o pangwasFilter.o pangwasIO.o pangwasAssoc.o

all: $(PROGRAMS)

clean:
	$(RM) *.o ~* $(PROGRAMS)

install: all
	install -d $(BINDIR)
	install $(PROGRAMS) $(BINDIR)

pangwas: $(CLASSES) $(OBJECTS)
	$(LINK.cpp) $^ $(LDLIBS) -o $@

.PHONY: all test clean install

