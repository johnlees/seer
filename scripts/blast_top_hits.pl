#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
use Text::CSV;

# Features
# BLAST and BLAT
# Exact/full query matches only
#

my $clean = 0;
my $tmp_sort = "tmp_blat.psl";

my $help_message = <<HELP;
Usage: ./blast_top_hits.pl --input blast_out --mode blast <options>

Reduces a blast or blat output file to the top hit for each query sequence
only

   Options

   Required
   --input              BLAST output file in -m 6/-outfmt 6 or BLAT .psl file
   --mode               Either blast or blat

   Optional
   --full_query        The entire query must be reported as a match to appear
                       in the output

   -h, --help          Shows this help.

HELP
#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($file_in, $mode, $full_query, $help);
GetOptions( "input|f=s" => \$file_in,
            "mode=s" => \$mode,
            "full_query" => \$full_query,
            "help|h" => \$help
		   ) or die($help_message);

# Basic error check on inputs
if (defined($help))
{
   print STDERR $help_message;
}
elsif (!defined($file_in) || !-e $file_in)
{
   print STDERR "Input file not found\n";
   print STDERR $help_message;
}
elsif (!defined($mode))
{
   print STDERR "Operation mode must be either blast or blat\n";
   print STDERR $help_message;
}
else
{
   my (@columns, $qseqid, $qstart, $qend, $length);
   if ($mode =~ /blat/i)
   {
      # psl format
      # match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts
      #         match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
      @columns = qw( match mismatch repmatch Ns Qgapcount Qgapbases Tgapcount Tgapbases strand Qname Qsize Qstart Qend Tname Tsize Tstart Tend blockcount blocksizes qstarts tstarts );

      $mode = "blat";
      $qseqid = "Qname";
      $qstart = "Qstart";
      $qend = "Qend";
      $length = "match";

      # Ensure sorted by query, then match score
      system("sort -k10,10 -k1,1nr $file_in > $tmp_sort");

      $file_in = $tmp_sort;
      $clean = 1;
   }
   elsif ($mode =~ /blast/i)
   {
      # blast outfmt 6
      #qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
      @columns = qw( qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore );

      $mode = "blast";
      $qseqid = "qseqid";
      $qstart = "qstart";
      $qend = "qend";
      $length = "length;"
   }
   else
   {
      die("Operation mode $mode not supported. Use either blast or blat\n");
   }

   my $csv = Text::CSV->new ( { binary => 1, sep_char => "\t", auto_diag => 1, eol => $/ } )  # should set binary attribute.
                 or die "Cannot use CSV: ".Text::CSV->error_diag ();
   my %rec;
   $csv->bind_columns (\@rec{@columns});

   open my $blast_in, "<:encoding(utf8)", "$file_in" or die "$file_in: $!";

   my $previous_query = "";
   while ($csv->getline($blast_in))
   {
      # Using speed advice from http://www.perlmonks.org/?node_id=937023
      if ($rec{$qseqid} ne $previous_query)
      {
         if (!$full_query || ((($mode eq "blast" && $rec{$qstart} == 1) || $rec{$qstart} == 0) && $rec{$qend} == $rec{$length}))
         {
            $csv->print(*STDOUT, [ @rec{@columns} ]);
         }
         $previous_query = $rec{$qseqid};
      }
   }

   if ($clean)
   {
      unlink $file_in;
   }
}

exit(0);

