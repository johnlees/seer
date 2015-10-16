PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

CXX=icpc
CXXFLAGS=-Wall -O3 -parallel -ipo -std=c++11
CPPFLAGS=-I$(PREFIX)/include -I../gzstream -I../dlib -D DLIB_NO_GUI_SUPPORT=1 -D DLIB_USE_BLAS=1 -D DLIB_USE_LAPACK=1
SEER_LDLIBS=-L../gzstream -L$(PREFIX)/lib -lhdf5 -lgzstream -lz -larmadillo -lboost_program_options -mkl
MAP_LDLIBS=-L$(PREFIX)/lib -lboost_program_options -mkl
COMBINE_LDLIBS=-L../gzstream -L$(PREFIX)/lib -lgzstream -lz -lboost_program_options
FILTER_LDLIBS=-L$(PREFIX)/lib -lboost_program_options -mkl

PROGRAMS=seer kmds map_back combineKmers filter_seer

CLASSES=sample.o kmer.o covar.o
COMMON_OBJECTS=$(CLASSES) seerCommon.o seerErr.o seerFilter.o seerIO.o seerChiFilter.o
SEER_OBJECTS=$(COMMON_OBJECTS) seerMain.o seerCmdLine.o seerContinuousAssoc.o seerBinaryAssoc.o logitFunction.o
KMDS_OBJECTS=$(COMMON_OBJECTS) kmdsMain.o kmdsStruct.o kmdsCmdLine.o
MAP_OBJECTS=fasta.o significant_kmer.o mapMain.o mapCmdLine.o
COMBINE_OBJECTS=combineInit.o combineCmdLine.o combineKmers.o
FILTER_OBJECTS=significant_kmer.o filter_seer.o filterCmdLine.o

all: $(PROGRAMS)

clean:
	$(RM) *.o ~* $(PROGRAMS)

install: all
	install -d $(BINDIR)
	install $(PROGRAMS) $(BINDIR)

seer: $(SEER_OBJECTS)
	$(LINK.cpp) $^ $(SEER_LDLIBS) -o $@

kmds: $(KMDS_OBJECTS)
	$(LINK.cpp) $^ $(SEER_LDLIBS) -o $@

map_back: $(MAP_OBJECTS)
	$(LINK.cpp) $^ $(MAP_LDLIBS) -o $@

combineKmers: $(COMBINE_OBJECTS)
	$(LINK.cpp) $^ $(COMBINE_LDLIBS) -o $@

filter_seer: $(FILTER_OBJECTS)
	$(LINK.cpp) $^ $(FILTER_LDLIBS) -o $@

.PHONY: all test clean install

