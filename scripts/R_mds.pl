#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;

my $usage_message = <<USAGE;
Usage: ./R_mds.pl -k <filtered_kmers.txt> -p <pheno_file>

Projects distance matrix into lower number of dimensions
Requires R and rhdf5 to be installed, see github wiki for more details

   Options

   Required
   -d, --dist       Distance matrix output by kmds
   -p, --pheno      Phenotype file used by seer

   Optional
   -o, --output     Output prefix (default 'all_structure')
   --pc             Number of dimensions to project into (default 3)
   -R               R executable (default 'which R')
   -h, --help       Displays this message

USAGE

my $tmp_file_name = "mds_tmp.Rscript";

#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#

#* gets input parameters
my ($distance_file, $pheno, $output_prefix, $dimensions, $R_exec, $help);
GetOptions ("dist|d=s"  => \$distance_file,
            "pheno|p=s" => \$pheno,
            "output|o=s" => \$output_prefix,
            "pc=i" => \$dimensions,
            "R=s" => \$R_exec,
            "help|h"     => \$help
		   ) or die($usage_message);

if (defined($help) || !defined($distance_file) || !defined($pheno))
{
   print $usage_message;
}
else
{
   # Set defaults
   if (!defined($dimensions))
   {
      $dimensions = 3;
   }
   if (!defined($R_exec))
   {
      $R_exec = `which R`;
      chomp $R_exec;
   }
   if (!defined($output_prefix))
   {
      $output_prefix = "all_structure";
   }

   # Check number of samples matches, write .samples if ok
   open(PHENO, $pheno) || die("Could not open pheno file $pheno: $!\n");

   my @samples;
   while (my $line_in = <PHENO>)
   {
      chomp $line_in;
      my @fields = split(/\t/, $line_in);

      push(@samples, $fields[1]);
   }

   my @samples_out = sort @samples;

   open(DIST, $distance_file) || die("Could not open distance file $distance_file: $!\n");
   my $dist_line = <DIST>;
   my @dist_cols = split(",", $dist_line);

   if(scalar(@dist_cols) != scalar(@samples_out))
   {
      die("Number of samples in $pheno does not match those in $distance_file\n");
   }

   my $dist_rows = 1;
   while (my $line_in = <DIST>)
   {
      $dist_rows++;
   }

   if($dist_rows != scalar(@dist_cols))
   {
      die("Distance matrix $distance_file is not square\n");
   }

   close DIST;

   open(SAMPLES, ">$output_prefix.samples") || die("Could not write to $output_prefix.samples\n");
   foreach my $sample (@samples_out)
   {
      print SAMPLES $sample . "\n";
   }
   close SAMPLES;

   #****************************************************************************************#
   #* Rscript                                                                              *#
   #****************************************************************************************#
   my $Rscript = <<RSCRIPT;
library(rhdf5)
distances <- read.csv("$distance_file", header=FALSE, stringsAsFactors=FALSE)
projection <- cmdscale(distances, k = $dimensions, eig = T)
covariates <- projection[["points"]]

for (i in 1:$dimensions)
{
   covariates[,i] <- covariates[,i]/max(abs(covariates[,i]))
}

h5createFile("$output_prefix")
h5write(covariates,"$output_prefix","dataset")
H5close()

pdf("scree_plot.pdf")
plot(projection[[2]][1:30], type='l', xlab="Dimensions", ylab="Eigenvalue", main="Scree plot")
dev.off()
RSCRIPT

   open(RSCRIPT, ">$tmp_file_name") || die("Could not open $tmp_file_name: $!\n");
   print RSCRIPT $Rscript;
   close(RSCRIPT);

   print STDERR "Running projection\n";
   my $R_command = "$R_exec CMD BATCH --no-save --no-restore $tmp_file_name";
   system($R_command);

   unlink($tmp_file_name, "$tmp_file_name.Rout");
}

exit(0);

