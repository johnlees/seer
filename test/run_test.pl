#!/usr/bin/perl -w

use strict;
use warnings;

use File::Temp qw/ :POSIX /;

my $exit_status = 0;
my $seer_location = "../src/";

sub do_test($$$$)
{
   my ($command, $file, $num, $name) = @_;

   my $fail = 0;

   my $outfile = tmpnam();
   my $errfile = tmpnam();
   system("$command > $outfile 2> $errfile");

   my $outdiff = `diff -q $outfile $file.o.txt`;
   my $errdiff = `diff -q $errfile $file.e.txt`;
   if ($outdiff ne "" || $errdiff ne "")
   {
      print STDERR "FAILED test $num, $name\n";
      print STDERR $outdiff . "\n";
      print STDERR $errdiff . "\n";
      $fail = 1;
   }
   else
   {
      print STDERR "PASSED test $num, $name\n";
   }

   unlink($outfile);

   return($fail);
}

$exit_status = do_test("$seer_location/seer -k example_kmers.gz -p subset.pheno", "test1", 1, "basic filters");
$exit_status= $exit_status || do_test("$seer_location/seer -k example_kmers.gz -p subset.pheno --pval 1 --chisq 1", "test2", 2, "binary phenotype assocation");
$exit_status = $exit_status || do_test("$seer_location/seer -k example_kmers.gz -p subset.pheno --pval 1 --chisq 1 --maf 0.1 --print_samples", "test3", 3, "print output");
$exit_status = $exit_status || do_test("$seer_location/seer -k example_kmers.gz -p subset.cont.pheno --pval 1 --chisq 1 --maf 0.1", "test4", 4, "continuous phenotype assocation");
$exit_status = $exit_status || do_test("$seer_location/seer -k example_kmers.gz -p example.pheno --struct all_structure_new --pval 1 --chisq 1", "test5", 5, "assocation with population structure");
$exit_status = $exit_status || do_test("$seer_location/seer -k example_kmers.gz -p subset.pheno --covar_file covariates.txt --covar_list 2q,3 --pval 1 --chisq 1", "test6", 6, "assocation with covariates");
$exit_status = $exit_status || do_test("$seer_location/filter_seer -k filter_in.txt --pos_beta", "test7", 7, "filter output");
$exit_status = $exit_status || do_test("$seer_location/map_back -k map_in.txt -r assembly_locations.txt --threads 1", "test8", 8, "map k-mers");

exit($exit_status);

