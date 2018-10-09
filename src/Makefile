# Default prefix. Can change this to system directories if needed
#PREFIX=${HOME}/software
#BINDIR=$(PREFIX)/bin

# Intel compiler - uncomment if you have icpc and mkl.
# gcc
#CXX=icpc
#CXXFLAGS=-Wall -O3 -parallel -ipo -std=c++11
#SEER_LDLIBS=-L../gzstream -L$(PREFIX)/lib -lhdf5 -lgzstream -lz -larmadillo -lboost_program_options -mkl
#MAP_LDLIBS=-L$(PREFIX)/lib -lboost_program_options -mkl
#COMBINE_LDLIBS=-L../gzstream -L$(PREFIX)/lib -lgzstream -lz -lboost_program_options
#FILTER_LDLIBS=-L$(PREFIX)/lib -lboost_program_options -mkl
# gcc
CXXFLAGS=-Wall -O3 -std=c++11
SEER_LDLIBS=-L../gzstream -L$(PREFIX)/lib -L/usr/local/hdf5/lib -lhdf5 -lgzstream -lz -larmadillo -lboost_program_options -llapack -lblas -lpthread
MAP_LDLIBS=-L$(PREFIX)/lib -lboost_program_options -lpthread
COMBINE_LDLIBS=-L../gzstream -L$(PREFIX)/lib -lgzstream -lz -lboost_program_options
FILTER_LDLIBS=-L$(PREFIX)/lib -lboost_program_options

CPPFLAGS=-I$(PREFIX)/include -I../gzstream -I../dlib -I/usr/local/hdf5/include -D DLIB_NO_GUI_SUPPORT=1 -D DLIB_USE_BLAS=1 -D DLIB_USE_LAPACK=1 -DARMA_USE_HDF5=1

# For static linking. To use type 'make static'
# Requires icpc and mkl
# Dynamically links libm, libz, libc
#COMMON_LDLIBS=-L../gzstream -L$(PREFIX)/lib -static-libstdc++ -static-libgcc -static-intel
#SEER_STATIC_LDLIBS=$(COMMON_LDLIBS) -Wl,-Bstatic -lhdf5 -lgzstream -larmadillo -lboost_program_options -Wl,-Bdynamic -lz -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group
#MAP_STATIC_LDLIBS=$(COMMON_LDLIBS) -Wl,-Bstatic -lboost_program_options -Wl,-Bdynamic
#COMBINE_STATIC_LDLIBS=$(COMMON_LDLIBS) -Wl,-Bstatic -lgzstream -lboost_program_options -Wl,-Bdynamic -lz
#FILTER_STATIC_LDLIBS=-$(COMMON_LDLIBS) -Wl,-Bstatic -lboost_program_options -Wl,-Bdynamic

# Full static linking with gcc
COMMON_LDLIBS=-L../gzstream -L$(PREFIX)/lib -static -static-libstdc++ -static-libgcc
SEER_STATIC_LDLIBS=$(COMMON_LDLIBS) -L../gzstream -L$(PREFIX)/lib -lhdf5 -lgzstream -lz -larmadillo -lboost_program_options -llapack -lblas -lpthread -lgfortran -lquadmath -lm
MAP_STATIC_LDLIBS=$(COMMON_LDLIBS) -L$(PREFIX)/lib -lboost_program_options -lpthread
COMBINE_STATIC_LDLIBS=$(COMMON_LDLIBS) -L../gzstream -L$(PREFIX)/lib -lgzstream -lz -lboost_program_options
FILTER_STATIC_LDLIBS=$(COMMON_LDLIBS) -L$(PREFIX)/lib -lboost_program_options


PROGRAMS=seer kmds map_back combineKmers filter_seer
STATIC_PROGRAMS=seer_static kmds_static map_back_static combineKmers_static filter_seer_static

CLASSES=sample.o significant_kmer.o kmer.o covar.o
COMMON_OBJECTS=$(CLASSES) seerCommon.o seerErr.o seerIO.o seerBasicFilter.o
SEER_OBJECTS=$(COMMON_OBJECTS) seerMain.o seerCmdLine.o seerStats.o seerContinuousAssoc.o seerBinaryAssoc.o logitFunction.o linearFunction.o
KMDS_OBJECTS=$(COMMON_OBJECTS) kmdsMain.o kmdsStruct.o kmdsCmdLine.o
MAP_OBJECTS=fasta.o significant_kmer.o mapMain.o mapCmdLine.o
COMBINE_OBJECTS=combineInit.o combineCmdLine.o combineKmers.o
FILTER_OBJECTS=significant_kmer.o filter_seer.o filterCmdLine.o

all: $(PROGRAMS)

static: $(STATIC_PROGRAMS)

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

seer_static: $(SEER_OBJECTS)
	$(LINK.cpp) $^ $(SEER_STATIC_LDLIBS) -o seer

kmds_static: $(KMDS_OBJECTS)
	$(LINK.cpp) $^ $(SEER_STATIC_LDLIBS) -o kmds

map_back_static: $(MAP_OBJECTS)
	$(LINK.cpp) $^ $(MAP_STATIC_LDLIBS) -o map_back

combineKmers_static: $(COMBINE_OBJECTS)
	$(LINK.cpp) $^ $(COMBINE_STATIC_LDLIBS) -o combineKmers

filter_seer_static: $(FILTER_OBJECTS)
	$(LINK.cpp) $^ $(FILTER_STATIC_LDLIBS) -o filter_seer


.PHONY: all static test clean install

