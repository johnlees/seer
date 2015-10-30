export PREFIX=${HOME}/software
export BINDIR=$(PREFIX)/bin

all:
	cd src && $(MAKE) all
	cd test && $(MAKE) all

clean:
	cd src && $(MAKE) clean
	cd test && $(MAKE) clean

install: all
	cd src && $(MAKE) install

test: all
	cd test && $(MAKE) test

.PHONY: all clean install test

