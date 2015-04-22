#!/usr/bin/env python

#
# PANPIPES.py
# Pangwas ANalysis PIPEline Script
#
#

import os,sys
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="A file containing a list of dsm output files")
parser.add_argument("pheno", help=".pheno file containing metadata")
parser.add_argument("-o", "--out_prefix", help="Prefix for output files",default="pangwas")
parser.add_argument("-t","--threads", help="Threads to use",type=int,default=4)
parser.add_argument("--LSF", help="Submit over LSF",action="store_true",default=False)
parser.add_argument("--pcs", help="Number of principal coordinates to use in population structure",type=int,default=3)
parser.add_argument("--subsample", help="Proportion of total kmers to use in mds calculation",type=float,default=0.001)
parser.add_argument("--maf", help="Minimum minor allele frequency",type=float,default=0.01)
parser.add_argument("--chi2",help="Chi^2 filter cutoff",type=float,default=10e-5)
parser.add_argument("--pval", help="pvalue cutoff",type=float,default=10e-8)
parser.add_argument("--assemble", help="assemble significant kmers, and perform association",action="store_true",default=False)
parser.add_argument("--reference", help="map kmers back to this reference file")
parser.add_argument("--drafts", help="file of annotated draft assemblies to map kmers back to")
args = parser.parse_args()


# use blat by default
# fall back to blast

