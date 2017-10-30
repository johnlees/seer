#export PREFIX=${HOME}/software
#export BINDIR=$(PREFIX)/bin

all: gzstream
	cd src && $(MAKE) all

gzstream:
	cd gzstream && $(MAKE) gzstream

clean:
	cd src && $(MAKE) clean
	cd test && $(MAKE) clean

install: all
	cd src && $(MAKE) install

test: all
	cd test && $(MAKE) test

.PHONY: all clean install test

