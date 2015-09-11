# seer
Sequence element enrichment analysis

Installation
==============
First clone the repository

    git clone https://github.com/johnlees/seer

Also include dlib, if necessary

    git clone --recursive https://github.com/johnlees/seer

Currently tested on Linux only, installation should proceed as

    cd src
    make
    make install

Dependencies
--------------
seer currently depends on

- boost <http://www.boost.org/>
- gzstream <http://www.cs.unc.edu/Research/compgeom/gzstream/>
- armadillo <http://arma.sourceforge.net/>
- dlib <http://www.dlib.net/>
- HDF5 <https://www.hdfgroup.org/HDF5/>

You will also require

- gcc >4.9 or equivalent
- gcc libstdc++ >4.9

Brief installation instructions

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

**gzstream**

Download and unpack. Change into directory

    make

**armadillo**

Download and unpack. Change into directory

    cmake .
    make
    make install

**dlib**

If not installed use the above git clone command to include with the
repository. Otherwise unpack header files to $(PREFIX)/include

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

Usage
=============
The wrapper script is crystal_ball.py, in the top level directory. Main programs are installed into bin.

Count your kmers
--------------
First you'll need to count the kmers which are present in your samples, which is best done from fasta files of assembled reads. If possible, we recommend the use of [distributed string mining](https://github.com/HIITMetagenomics/dsm-framework) as this is memory efficient and selects kmer size based on entropy, rather than being a fixed arbitrary word length.

If this is not possible, use [dsk](http://minia.genouest.org/dsk/) to count your kmers, followed by the combineKmers program we provide:

    dsk -file sample1_contigs.fa -abundance-min 1 -out sample1_dsk
    dsk2ascii -file sample1_dsk.h5 -out sample1_dsk.txt
    ./combineKmers -r 21mer_samples.txt -o all_21mers_new --min_samples 2

The final step is fairly memory intensive. Taking the union of 21mers in 3000 samples took around 25Gb of RAM and about 90 minutes, however this step only needs to be completed once.

If using dsk, we've found using kmers between 21 and 41 seems to work well.

The output of both these methods can be read directly by kmds and seer.

kmds
--------------
Next, filter kmers and create a matrix representing population structure with kmds

    ./kmds -k dsm_input.txt.gz --pheno metadata.pheno -o filtered

To spread this process out, run the following command on each dsm file

    ./kmds -k dsm_input.txt.gz --pheno metadata.pheno --no_mds --size 10000

Then use

    ./kmds --mds_concat subsampled_matrices.txt -o all_structure --threads 16

seer
--------------
Type ./seer with no options to get brief usage, or for a full option listing

    ./seer -h

Basic usage is as follows

    ./seer -k dsm_input.txt.gz --pheno metadata.pheno > significant_kmers.txt

To use the kmds output, and increase execution speed

    ./seer -k filtered.gz --pheno metadata.pheno --struct filtered.dsm --no_filtering --threads 4

The same filtering options as in kmds are available, as well a filtering on the adjusted p-value

Output columns are sequence of element, frequency in sample set,
uncorrected p-value, corrected p-value, odds-ratio of effect size and
optionally samples the sequence element appears in

Interpreting the output
--------------
Run the post analysis program filter_seer

    ./filter_seer -k significant_kmers.txt --pos_beta --maf 0.05 --sort pval > seer.filtered.txt

Run map_back to find the position of significant kmers, in the assemblies they are from

    ./map_back -k seer.filtered.k31.txt -r references.txt --threads 4

Run blat or blast on the significant kmers

    cut -f 1 seer.filtered.k31.txt | awk '{print ">"NR"\n"$1}' > seer.sorted.k31.fa
    blat -noHead -stepSize=2 -minScore=15 reference.fa seer.sorted.k31.fa seer.k31.blat.psl
    blast -task blastn-short -subject reference.fa -query seer.sorted.k31.fa -outfmt 6 seer.k31.blast.txt

Reduce this to the top hits

    ./scripts/blast_top_hits.pl --input seer.k31.blat.psl --mode blat --full_query > top_hits.txt

If you have [bedops](http://bedops.readthedocs.org/en/latest/index.html) and [bedtools](http://bedtools.readthedocs.org/en/latest/) you can annotate the position of these hits, for example

    psl2bed < seer.k31.blat.psl > seer.k31.blat.bed
    gff2bed < reference.gff > reference.bed
    bedtools intersect -a seer.k31.blat.bed -b reference.bed -wb
    bedtools intersect -a seer.k31.blat.bed -b reference.bed -c | sort -n -r -t $'\t' -k11,11


Troubleshooting
=============

General
-------------
**When I run the programs I get an error about a missing GLIBCXX**

    seer: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.17' not found (required by seer)
    seer: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by seer)

Make sure you've got the gcc libstdc++ v4.9 or higher available. The
default on your OS may be a much lower version than this

**Analysis is slow**

Most steps can be effectively parallelised in two ways. Increase the threading of the step with --threads, and/or split the input files and run each command independently on each input file.

    split -n l/15 all_31mers 31mers

will split the input file all_31mers into 15 pieces which can be analysed independently. Use 'cat' to combine the results.

kmds
-------------
**This takes ages to run**

First, split up your input files into say 16 pieces, then subsample each
one in its own process using --no_mds and --size.

Then list these output files in a text file, and input using
--mds_concat. Using 16 threads this took about 30 hour and 53Gb on our test
system of 3 000 samples, using 1% of the kmers.

Choose the number of kmers to subsample carefully, lower and this will run
faster and require less memory. A size of 1% of the number of kmers is
appropriate, but as low as 0.1% should work.

seer
-------------
**The p-values are slightly different from what I expect**

For very small p-values (when W > 5) an upper bound is calculated rather
than the exact p-value for robustness and computational speed. This
bound is accurate to +/- W^(-2) %

If you need the exact value, checkout the erfc branch which uses
arbitrary precision floats and doesn't rely on this bound.

**I get p-values of 0**

p-value is too small to represent as a float (<10E-308). At this point
it is a bit meaningless to assign a value, but essentially the data
exactly fits the regression.

Biologically, either the effect size is enormous, or more likely the
phenotype is monophyletic (see below).

**I get the error message: 'inv_sympd(): matrix appears to be singular'**

The optimiser may not have converged, check previous error messages and
see the note about this below.

If it has, it's likely the p-value is too small to calculate. You can
approximate this p-value as 0 if you wish (the program default). If
you're getting lots of messages like this your phenotype is probably
perfectly associated with a clade, or small number of clades.

Check whether this is the case using a tool such as [phylocanvas](http://phylocanvas.org/).
If your phenotype is monophyletic then you won't have good resolution on
which sequence elements may be related to it. See the paper for more
explanation on this

**I get the error message: 'kmer convergence error: The objective function generated non-finite outputs'**

The mds structure is probably poorly scaled compared to the kmers. kmers
presence and absence is coded as 1 and 0 respectively.

kmds produced a rectangular matrix of mds values, where each row is
a sample, and each column is a decreasing dimension. Each dimension is
scaled to values in the range [-1,1].

These matrices are stored in hdf5 format, so you can use these tools to
inspect them, and if necessary rescale them. Most convenient is the
R package [rhdf5](http://www.bioconductor.org/packages/release/bioc/html/rhdf5.html).

