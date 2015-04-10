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

You will also require

- gcc >4.8 or equivalent

Brief installion instructions

**boost**

Best installed with your distribution's package manager, and you should use the c++11 version if possible.

For a manual installation, see <http://www.boost.org/doc/libs/1_57_0/more/getting_started/unix-variants.html> for details on how to use ./b2 to install. I reccommend that you create a user-config.jam file in the boost root which modifies the gcc compilation:

   using gcc:
      : std11
      : g++
      : <cxxflags>-std=c++11

Then run

   ./bootstrap.sh
   ./b2 install toolset=gcc-std11

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
Type ./pangwas with no options to get brief usage, or for a full option listing

    ./pangwas -h

Basic usage is as follows

    ./pangwas -k dsm_input.txt --pheno metadata.pheno > significant_kmers.txt

