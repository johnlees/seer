# pangwas
Implementation of a pan-genome wide association study by using kmers

Installation
==============
First clone the repository
   git clone https://github.com/johnlees/pangwas

Also include mlpack, if necessary
   git clone --recursive https://github.com/johnlees/pangwas

Currently tested on Linux only, installation should proceed as
   cd src
   make
   make install

Dependencies
--------------
pangwas currently depends on
- boost <http://www.boost.org/>
- armadillo <http://arma.sourceforge.net/>
- mlpack <http://www.mlpack.org/>

Brief installion instructions

**boost**
TODO
Should use -std=c++11

**armadillo**
Download and unpack. Change into directory
   cmake .
   make
   make install

**mlpack**
Download and unpack. Change into directory
   mkdir build && cd build
   cmake ../
   make logistic_regression
   make install logistic_regression

Usage
=============
TODO
