#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $usage_message = <<USAGE;
Usage: ./hits_to_matrix.pl -k <filtered_kmers.txt> -p <pheno_file>

Reformats output for easier reading as tsv. Prints to stdout

   Options

   Required
   -k, --kmers      Filtered significant kmers, as output by seer + filter_seer
   -p, --pheno       Phenotype file used by seer

   Optional
   -h, --help       Displays this message

USAGE

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($kmer_file, $pheno, $help);
GetOptions ("kmers|k=s"  => \$kmer_file,
            "pheno|p=s" => \$pheno,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help) || !defined($kmer_file) || !defined($pheno))
{
   print $usage_message;
}
else
{
   # Read possible samples
   open(PHENO, $pheno) || die("Could not open pheno file $pheno\n");

   my @samples;
   while (my $pheno_line = <PHENO>)
   {
      chomp $pheno_line;

      my ($sample, $junk, $pheno) = split("\t", $pheno_line);

      push(@samples, $sample);
   }

   # Keep everything sorted (as depends on kmer counting program used)
   @samples = sort @samples;
   close PHENO;

   # Print header row to stdout
   open(KMERS, $kmer_file) || die("Could not open kmer file $kmer_file\n");
   print join("\t", "sequence", "maf", "unadj_p_val", "p_val", "beta", "se", "comments", @samples) . "\n";

   while (my $kmer_line = <KMERS>)
   {
      chomp $kmer_line;

      my ($sequence, $maf, $unadj, $adj, $beta, $se, $comment, @kmer_samples) = split("\t", $kmer_line);
      @kmer_samples = sort @kmer_samples;

      print join("\t", $sequence, $maf, $unadj, $adj, $beta, $se, $comment) . "\t";

      my $next_sample = shift(@kmer_samples);
      foreach my $sample (@samples)
      {
         print "\t";
         if ($sample eq $next_sample)
         {
            print $sample;
            if (scalar(@kmer_samples) != 0)
            {
               $next_sample = shift(@kmer_samples);
            }
         }
      }
      print "\n";

   }
   close KMERS;
}
