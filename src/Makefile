PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

CXX=~/homebrew/bin/g++-4.9
CXXFLAGS=-Wall -g -std=c++11
#CXXFLAGS=-Wall -g -std=c++11 -stdlib=libc++
CPPFLAGS=-I../armadillo-4.650.3/include/ -I../mlpack-1.0.12/build/include/ -I$(PREFIX)/include -I/usr/include/libxml2
LDLIBS=-L../mlpack-1.0.12/build/lib -L../armadillo-4.650.3 -L$(PREFIX)/lib -lm -framework Accelerate -lmlpack -larmadillo -lboost_program_options

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

