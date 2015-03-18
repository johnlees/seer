PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

BOOSTDIR=${HOME}/homebrew

CXXFLAGS=-Wall -g -std=c++11 -stdlib=libc++
CPPFLAGS=-I$(PREFIX)/include -I$(BOOSTDIR)/include -I/usr/include/libxml2
LDLIBS=-L$(PREFIX)/lib -L$(BOOSTDIR)/lib -lm -lboost_program_options -larmadillo -lmlpack

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

