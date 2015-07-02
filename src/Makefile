PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

CXX=icpc
CXXFLAGS=-Wall -O3 -parallel -ipo -std=c++11
CPPFLAGS=-I$(PREFIX)/include -I../gzstream -I../dlib -D DLIB_NO_GUI_SUPPORT=1 -D DLIB_USE_BLAS=1 -D DLIB_USE_LAPACK=1
LDLIBS=-L../gzstream -L$(PREFIX)/lib -lhdf5 -lgzstream -lz -larmadillo -lboost_program_options -mkl

PROGRAMS=pangwas pangloss map_back

CLASSES=sample.o kmer.o
COMMON_OBJECTS=$(CLASSES) panCommon.o panErr.o panFilter.o panIO.o panChiFilter.o
PANGWAS_OBJECTS=$(COMMON_OBJECTS) pangwasMain.o pangwasCmdLine.o pangwasContinuousAssoc.o pangwasBinaryAssoc.o logitFunction.o
PANGLOSS_OBJECTS=$(COMMON_OBJECTS) panglossMain.o panglossStruct.o panglossCmdLine.o
MAP_OBJECTS=fasta.o significant_kmer.o mapMain.o mapCmdLine.o

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

map_back: $(MAP_OBJECTS)
	$(LINK.cpp) $^ $(LDLIBS) -o $@

.PHONY: all test clean install

