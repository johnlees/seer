PREFIX=${HOME}/software
BINDIR=$(PREFIX)/bin

BOOSTDIR=${HOME}/homebrew

CXXFLAGS=-Wall -g
CPPFLAGS=-I$(BOOSTDIR)/include
LDLIBS=-L$(BOOSTDIR)/lib -lm -lboost_program_options

PROGRAMS=pangwas

all: $(PROGRAMS)

clean:
	$(RM) *.o ~* $(PROGRAMS)

install: all
	install -d $(BINDIR)
	install $(PROGRAMS) $(BINDIR)

pangwas: pangwasMain.o pangwasFilter.o pangwasIO.o pangwasAssoc.o
	$(LINK.cpp) $^ $(LDLIBS) -o $@

.PHONY: all test clean install

