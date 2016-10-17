#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $usage_message = <<USAGE;
Usage: ./mapping_to_phandango.pl -b <blast_file> -a <fastlmm_kmers> > phandango.plot

Creates a plot file for visualisation in phandango

   Options

   Required
   -b, --blast      blast results
   -a, --assoc      fastlmm results

   -h, --help       Displays this message

USAGE

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($assoc_file, $blast_file, $help);
GetOptions ("assoc|a=s"  => \$assoc_file,
            "blast|b=s" => \$blast_file,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help) || !defined($assoc_file) || !defined($blast_file))
{
   print $usage_message;
}
else
{
   open(ASSOC, $assoc_file) || die("Could not open assoc file $assoc_file: $!\n");

   my @pvals;
   while (my $assoc_line = <ASSOC>)
   {
      chomp $assoc_line;

      my (@assoc_fields) = split("\t", $assoc_line);
      push(@pvals, -log($assoc_fields[5])/log(10));
   }

   close ASSOC;

   open(BLAST, $blast_file) || die("Could not open blast file $blast_file: $!\n");

   my %points;
   while (my $blast_line = <BLAST>)
   {
      chomp $blast_line;

      my (@blast_fields) = split("\t", $blast_line);

      # Get position, check for reverse complement
      my $start = $blast_fields[8];
      my $end = $blast_fields[9];

      if ($start > $end)
      {
         my $tmp = $start;
         $start = $end;
         $end = $tmp;
      }

      $points{$blast_fields[0]}{pval} = $pvals[$blast_fields[0]-1];
      $points{$blast_fields[0]}{start} = $start;
      $points{$blast_fields[0]}{pos} = "$start..$end";
   }

   # Sort and print the output
   my @keys = sort { $points{$a}{start} <=> $points{$b}{start} } keys(%points);
   foreach my $kmer (@keys)
   {
      print join("\t", "26", ".", $points{$kmer}{pos}, $points{$kmer}{pval}, "0") . "\n";
   }
}

exit(0);


