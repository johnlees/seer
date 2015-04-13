#!/usr/bin/perl -w

use strict;
use warnings;

my $sample_file = $ARGV[0];
my $sites = $ARGV[1];

my (@samples, %phenotypes);
open(PHENO, $sample_file) || die("Could not open $sample_file\n");
while (my $pheno_line = <PHENO>)
{
   chomp $pheno_line;
   my ($FID, $IID, $pheno) = split("\t", $pheno_line);

   push(@samples, "$IID:1");
   $phenotypes{"$IID:1"} = $pheno;
}
close PHENO;

# Sites increase in minor allele count (mac)
for (my $site_idx = 1; $site_idx <= $sites; $site_idx++)
{
   my @sampled_samples;

   my $mac = int(scalar(@samples) * ($site_idx/$sites));

   # Reservoir sampler over all samples for each site
   for (my $sample_idx = 0; $sample_idx < scalar(@samples); $sample_idx++)
   {
      if (scalar(@sampled_samples) < $mac)
      {
         push(@sampled_samples, $samples[$sample_idx]);
      }
      else
      {
         my $r = int(rand($sample_idx));
         if ($r < $mac)
         {
            $sampled_samples[$r] = $samples[$sample_idx];
         }
      }
   }

   print join(" ", $site_idx, "5.033430 0.161246 100 0 100 0.151841 100 |", @sampled_samples) . "\n";
}

exit(0);

