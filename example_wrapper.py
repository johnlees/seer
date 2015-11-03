#!/usr/bin/env python

#
# This is an example of a wrapper script for the whole seer analysis pipeline,
# based on having a cluster with LSF
#

# imports
import os,sys
from __future__ import print_function
import re
import argparse
import subprocess

# globals
subsampled_list = "kmds_tmp_files.txt"
job_num = re.compile('^Job <(\d+)>')

# subroutines

# Checks, for a list of job numbers, whether they are all done (or failed)
def check_done(jobs):
    done = 1
    for job in jobs:
        bsub_ret = subprocess.check_call("bjobs -a -noheader -o \"stat exit_code delimiter=','\" " + str(job), shell=True)

        status = bsub_ret.split(",")
        if bsub_ret[1] == "RUN" | bsub_ret[1] == "PEND" | bsub_ret[1] == "WAIT":
            done = 0
            break

    return done

# Command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("infile", help="A file containing a list of dsm output files")
parser.add_argument("pheno", help=".pheno file containing metadata")
parser.add_argument("-o", "--out_prefix", help="Prefix for output files",default="seer")
parser.add_argument("-t", "--threads", help="Threads to use",type=int,default=4)
parser.add_argument("--LSF", help="Submit over LSF",action="store_true",default=False)
parser.add_argument("--pcs", help="Number of principal coordinates to use in population structure",type=int,default=3)
parser.add_argument("--subsample", help="Proportion of total kmers to use in mds calculation",type=float,default=0.001)
parser.add_argument("--maf", help="Minimum minor allele frequency",type=float,default=0.01)
parser.add_argument("--chisq",help="Chi^2 filter cutoff",type=float,default=10e-5)
parser.add_argument("--pval", help="pvalue cutoff",type=float,default=10e-8)
parser.add_argument("--assemble", help="assemble significant kmers, and perform association",action="store_true",default=False)
parser.add_argument("--reference", help="map kmers back to this reference file")
parser.add_argument("--drafts", help="file of annotated draft assemblies to map kmers back to")
args = parser.parse_args()

# Read in dsm files
if not (os.path.isfile(args.infile) & os.path.isfile(args.pheno)):
    raise Exception("Mandatory input files do not exist")

with open(args.infile,'r') as f:
    dsm_files = f.readlines()

dsm_files = [x.strip('\n') for x in dsm_files]

# Run kmds --no_mds on each file in parallel
# Collect output
print("Filtering " + str(len(dsm_files)) + " files\n")
i = 0
subsampled_output = []
jobs = []

for dsm in dsm_files:
    i += 1

    length = args.subsample * int(subprocess.check_output("gzip -d -c " + dsm + " | wc -l", shell=True))
    try:
        kmds_command = ""
        if args.LSF:
            kmds_command = "bsub -o kmds.step1.%J." + str(i) + ".o -e kmds.step1.%J." + str(i) + ".e "

        kmds_command += "kmds -k " + str(dsm) + " -p " + str(args.pheno) + " -o kmds.step1." + str(i) + " --no_mds --maf " + str(args.maf) + " --chisq " + str(args.chi2) + " --size " + str(length)

        print(kmds_command)

        if args.LSF:
            job_return = subprocess.check_output(kmds_command, shell=True)
            m = job_num.match(job_return)

            jobs.append(m.group())
        else:
            retcode = subprocess.call(kmds_command, shell=True)
            if retcode < 0:
                print("kmds step 1 file " + str(i) + " failed with ", -retcode, file=sys.stderr)

        subsampled_output.append("kmds.step1." + str(i))
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

# Check all jobs have finished
if args.LSF:
    while not (check_done(jobs)):
        os.sleep(30)
    jobs = []

write_list = open(subsampled_list,'w')
for subsample in subsampled_output:
    write_list.write(subsample + "\n")

write_list.close()

# Run kmds --mds_concat on output
print("Calculating MDS components\n")
try:
    kmds_command = ""
    if args.LSF:
        kmds_command = "bsub -o kmds.step2.%J.o -e kmds.step2.%J..e -n" + str(args.threads) + " -R \"span[hosts=1]\" -R \"select[mem>4000] rusage[mem=4000]\" -M4000 "

    kmds_command += "kmds --mds_concat " + str(subsampled_list) + " -o all_structure --threads " + str(args.threads) + " --pc " + str(args.pcs)

    print(kmds_command)

    if args.LSF:
        job_return = subprocess.check_output(kmds_command, shell=True)
        m = job_num.match(job_return)

        jobs.append(m.group())
    else:
        retcode = subprocess.call(kmds_command, shell=True)
        if retcode < 0:
            print("kmds step 2 file failed with ", -retcode, file=sys.stderr)
except OSError as e:
    print("Execution failed:", e, file=sys.stderr)

if args.LSF:
    while not (check_done(jobs)):
        os.sleep(30)
    jobs = []

# Run seer on each file in parallel
print("Association on " + str(len(dsm_files)) + " files\n")
i = 0
seer_output = []
for dsm in dsm_files:
    i += 1

    try:
        seer_command = ""
        if args.LSF:
            seer_command = "bsub -o kmds.step1.%J." + str(i) + ".o -e kmds.step1.%J." + str(i) + ".e -n" + str(args.threads) + " -R \"span[hosts=1]\""

        seer_command += "'seer -k " + str(dsm) + " -p " + str(args.pheno) + " --no_filtering --pval " + str(args.pval) + " --struct all_structure.dsm --threads " + str(args.threads) + " --print_samples > seer." + str(i) + ".result'"

        print(seer_command)

        if args.LSF:
            job_return = subprocess.check_output(seer_command, shell=True)
            m = job_num.match(job_return)

            jobs.append(m.group())
        else:
            retcode = subprocess.call(seer_command, shell=True)
            if retcode < 0:
                print("seer file " + str(i) + " failed with ", -retcode, file=sys.stderr)

        seer_output.append("kmds." + str(i))
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)

if args.LSF:
    while not (check_done(jobs)):
        os.sleep(30)
    jobs = []

# TODO post-processing:
subprocess.check_call("cat seer.*.result > seer.result", shell=True)
for i in range(1, len(dsm_files)):
    subprocess.check_call("rm seer." + str(i) + ".result", shell=True)

# map back to references
subprocess.call("map_back -k seer.result -r " + args.drafts + " > kmer_draft_locations.txt")

# assembly of significant kmers
if (args.assemble):
    subprocess.call("VelvetOptimiser.pl --s 11 --e 71 --t 1 -f '-short -fasta ../../resultErytrhomycin.fa' --m -1 --k max --c tbp")

# use blat by default
# fall back to blast

print("All pans piped\nAssocation results written to seer.result")

