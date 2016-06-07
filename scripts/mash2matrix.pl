#!/usr/bin/perl -w

use strict;
use warnings;

my $usage_message = <<USAGE;
Usage: ./mash2matrix.pl <matrix_in> > all_distances.csv

Uses mash (http://mash.readthedocs.io/) to create all pairwise distances
between samples.

NB: the input to mash should be the same order as the .pheno file
generate <assemblies.fa> with:
cut -f 1 metadata.pheno | tr '\n' ' '

Then run
mash sketch -o reference <assemblies.fa>
mash dist reference.msh reference.msh > mash_distances.txt
./mash2matrix.pl mash_distances.txt > all_distances.csv

USAGE

my $dist_file = $ARGV[0];
if (!-e $dist_file)
{
   print STDERR $usage_message;
}
else
{
   open(INPUT, $dist_file) || die("Could not open $dist_file\n");

   my %row_idx;
   my @mat_out;
   my $i = 0;
   while (my $line_in = <INPUT>)
   {
      chomp $line_in;
      my ($el1, $el2, $dist, @junk) = split("\t", $line_in);

      if (!defined($row_idx{$el1}))
      {
         $row_idx{$el1} = $i++;
      }
      if (!defined($row_idx{$el2}))
      {
         $row_idx{$el2} = $i++;
      }

      $mat_out[$row_idx{$el1}][$row_idx{$el2}] = $dist;
   }

   for (my $i = 0; $i < scalar(@{$mat_out[0]}); $i++)
   {
      print join(",", @{$mat_out[$i]}) . "\n";
   }

}

exit(0);

