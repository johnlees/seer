#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $usage_message = <<USAGE;
Usage: ./hits_to_fastq.pl -k <significant_kmers.txt> -b <bonferroni_correction>

Creates a fastq for mapping from significant kmers

   Options

   Required
   -k, --kmers      Significant kmers, as output by seer
   -p, --bonf       Bonferroni correction to use (try 10e-8)

   Optional
   -h, --help       Displays this message

USAGE

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($kmer_file, $bonf, $help);
GetOptions ("kmers|k=s"  => \$kmer_file,
            "bonf|b=s" => \$bonf,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help) || !defined($kmer_file) || !defined($bonf))
{
   print $usage_message;
}
else
{
   open(KMERS, $kmer_file) || die("Could not open kmer file: $kmer_file\n");
   my $header = <KMERS>;

   my $i = 1;
   while (my $line_in = <KMERS>)
   {
      chomp $line_in;
      my ($kmer, $maf, $p_unadj, $p_adj, @junk) = split("\t", $line_in);

      my $corrected_p = $p_adj / $bonf;
      if ($corrected_p > 1)
      {
         $corrected_p = 1;
      }

      my $ascii_val = int(-10*log($corrected_p)/log(10)) + 33;
      if ($ascii_val < 33)
      {
         $ascii_val = 33;
      }
      elsif ($ascii_val > 126)
      {
         $ascii_val = 126;
      }

      my $qual = chr($ascii_val) x length($kmer);

      print join("\n", "\@$i", $kmer, "+", $qual) . "\n";

      $i++;
   }

   close KMERS
}

exit(0);

