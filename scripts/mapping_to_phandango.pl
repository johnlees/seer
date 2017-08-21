#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my %rev_map = (
   A => "T",
   T => "A",
   C => "G",
   G => "C",
);

sub rev_comp($)
{
   my ($seq) = @_;
   my @chars = split(//, $seq);

   my $rev;
   for (my $i = scalar(@chars) - 1; $i < 0; $i--)
   {
      my $rev .= $rev_map{$chars[$i]};
   }

   return $rev;
}

my $usage_message = <<USAGE;
Usage: ./mapping_to_phandango.pl -s <sam_file> -k <seer_kmers> > phandango.plot

Creates a plot file for visualisation in phandango

   Options

   Required
   -s, --sam        Sam file of kmers mapped to reference (unsorted)
   -k, --kmers      Kmers that were mapped

   Optional
   -m, --map        Minimum mapping quality
   --print_kmer     Also print k-mer sequence in output

   -h, --help       Displays this message

USAGE

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($kmer_file, $print_kmer, $sam_file, $min_qual, $help);
GetOptions ("kmers|k=s"  => \$kmer_file,
            "print_kmer" => \$print_kmer,
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
   my @header_fields = split("\t", $header);
   my $lrt = 0;
   if ($header_fields[4] eq "lrt_p_val")
   {
      $lrt = 1;
   }

   my %points;
   while (my $sam_line = <SAM>)
   {
      chomp $sam_line;

      if ($sam_line !~ /^@/)
      {
         my ($kmer_id, $flag, $chrom, $pos, $map_qual, $cigar, $mate_chrom, $mate_pos, $insert_size, $sequence, $quality, @optional) = split("\t", $sam_line);

         # Get position
         my $start = $pos;
         my $end = $pos + length($sequence) - 1;

         my $kmer = <KMERS>;
         chomp $kmer;

         # Don't report unmapped or low quality mapped kmers
         if ($map_qual > $min_qual && !($flag & 4))
         {
            my ($kmer, $maf, $unadj, $adj, @other) = split("\t", $kmer);
            if ($lrt)
            {
               $adj = $other[0];
            }

            if ($kmer ne $sequence && $kmer ne rev_comp($sequence))
            {
               print STDERR "WARNING: line $kmer_id sequence mismatch. Are input files in same order?\n";
            }

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
      if (defined($print_kmer))
      {
         print join("\t", "26", $kmer, $points{$kmer}{pos}, $points{$kmer}{pval}, "0") . "\n";
      }
      else
      {
         print join("\t", "26", ".", $points{$kmer}{pos}, $points{$kmer}{pval}, "0") . "\n";
      }
   }
}

exit(0);

