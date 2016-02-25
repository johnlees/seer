#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $usage_message = <<USAGE;
Usage: ./mapping_to_phandango.pl -s <sam_file> -k <seer_kmers> > phandango.plot

Creates a plot file for visualisation in phandango

   Options

   Required
   -s, --sam        Sam file of kmers mapped to reference (unsorted)
   -k, --kmers      Kmers that were mapped

   Optional
   -m, --map        Minimum mapping quality

   -h, --help       Displays this message

USAGE

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($kmer_file, $sam_file, $min_qual, $help);
GetOptions ("kmers|k=s"  => \$kmer_file,
            "sam|s=s" => \$sam_file,
            "map|m=s" => \$min_qual,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help) || !defined($kmer_file) || !defined($sam_file))
{
   print $usage_message;
}
else
{
   if (!defined($min_qual))
   {
      $min_qual = 0;
   }

   open(SAM, $sam_file) || die("Could not open sam file $sam_file: $!\n");
   open(KMERS, $kmer_file) || die("Could not open kmer file $kmer_file: $!\n");

   my $header = <KMERS>;

   my %points;
   while (my $sam_line = <SAM>)
   {
      chomp $sam_line;

      if ($sam_line !~ /^@/)
      {
         my ($kmer_id, $flag, $chrom, $pos, $map_qual, $cigar, $mate_chrom, $mate_pos, $insert_size, $sequence, $quality, @optional) = split("\t", $sam_line);

         # Get position, check for reverse complement
         my ($end, $start);
         if ($flag & 16)
         {
            $end = $pos;
            $start = $pos - length($sequence) + 1;
         }
         else
         {
            $start = $pos;
            $end = $pos + length($sequence) - 1;
         }

         my $kmer = <KMERS>;
         chomp $kmer;

         # Don't report unmapped or low quality mapped kmers
         if ($map_qual > $min_qual && !($flag & 4))
         {
            my ($kmer, $maf, $unadj, $adj, @other) = split("\t", $kmer);

            my $log_p = 386; # Exponent limit of a double
            if ($adj > 0)
            {
               $log_p = -log($adj)/log(10);
            }

            $points{$kmer}{pval} = $log_p;
            $points{$kmer}{start} = $start;
            $points{$kmer}{pos} = "$start..$end";
         }
      }
   }

   # Sort and print the output
   my @keys = sort { $points{$a}{start} <=> $points{$b}{start} } keys(%points);
   foreach my $kmer (@keys)
   {
      print join("\t", "26", ".", $points{$kmer}{pos}, $points{$kmer}{pval}, "0") . "\n";
   }
}

exit(0);

