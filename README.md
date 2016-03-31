# seer
Sequence element enrichment analysis. This document contains
installation instuctions. Usage can be found on the [wiki](https://github.com/johnlees/seer/wiki/Usage), and more information in the [paper](http://biorxiv.org/content/early/2016/03/02/038463).

Installation
==============
First clone the repository

    git clone --recursive https://github.com/johnlees/seer

If you already have dlib:

    git clone https://github.com/johnlees/seer

Currently tested on Linux only, installation should proceed as

    make
    make install

Dependencies
--------------
seer currently depends on

- gzstream <http://www.cs.unc.edu/Research/compgeom/gzstream/>
- armadillo <http://arma.sourceforge.net/>
- boost <http://www.boost.org/>
- dlib <http://www.dlib.net/>
- HDF5 <https://www.hdfgroup.org/HDF5/>

You will also require

- gcc >4.9 or equivalent
- gcc libstdc++ >4.9

You probably already have boost, HDF5 and dlib (as long as you did clone --recursive).

###Brief installation instructions

**gzstream**

Download and unpack to a folder gzstream in the root of the repository. Change into the directory and type

    make

**HDF5**

Best installed with your distribution's package manager. Otherwise use
a binary from <https://www.hdfgroup.org/HDF5/release/obtain5.html>, or
if you wish to compile from source

    gunzip < hdf5-X.Y.Z.tar.gz | tar xf -
    cd hdf5-X.Y.Z
    ./configure --prefix=/usr/local/hdf5 <more configure_flags>
    make
    make check
    make install
    make check-install

**armadillo**

Make sure HDF5 is installed first.

Download and unpack. Change into directory and type

    cmake -DARMA_USE_HDF5=1 .
    make
    make install

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


**dlib**

If not installed use the above git clone command to include with the
repository. Otherwise unpack header files to $(PREFIX)/include

**installation**

Currently tested on Linux only, installation should proceed as

    make
    make install

You may need to explicitly set the current GCC compiler, which you can
do by running

    make CXX=g++-4.9


Usage, interpretation of results, and troubleshooting
=============
See the [wiki](https://github.com/johnlees/seer/wiki/Usage)
