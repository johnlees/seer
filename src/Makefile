PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

CXX=icpc
CXXFLAGS=-Wall -g -O3 -std=c++11
CPPFLAGS=-I../armadillo-4.650.3/include -I../mlpack/build/include -I$(PREFIX)/include -I/usr/include/libxml2 -I../gzstream
LDLIBS=-L../mlpack/build/lib -L../armadillo-4.650.4 -L../gzstream -L$(PREFIX)/lib -lgzstream -lz -lmlpack -larmadillo -lboost_program_options -mkl

#OS X compile options:
#CXX=clang++
#CXXFLAGS=-Wall -g -std=c++11 -stdlib=libc++
#LDLIBS=-L../mlpack-1.0.12/build/lib -L../armadillo-4.650.3 -L$(PREFIX)/lib -lm -framework Accelerate -lmlpack -larmadillo -lboost_program_options

PROGRAMS=pangwas pangloss

CLASSES=sample.o kmer.o
COMMON_OBJECTS=$(CLASSES) panCommon.o panErr.o panFilter.o panIO.o panChiFilter.o
PANGWAS_OBJECTS=$(COMMON_OBJECTS) pangwasMain.o pangwasCmdLine.o pangwasAssoc.o
PANGLOSS_OBJECTS=$(COMMON_OBJECTS) panglossMain.o panglossStruct.o panglossCmdLine.o

all: $(PROGRAMS)

clean:
	$(RM) *.o ~* $(PROGRAMS)

install: all
	install -d $(BINDIR)
	install $(PROGRAMS) $(BINDIR)

pangwas: $(PANGWAS_OBJECTS)
	$(LINK.cpp) $^ $(LDLIBS) -o $@

pangloss: $(PANGLOSS_OBJECTS)
	$(LINK.cpp) $^ $(LDLIBS) -o $@

.PHONY: all test clean install

