/*
 * File: sample.cpp
 *
 * Helper functions for the sample class
 *
 */

#include "sample.hpp"

const int default_pheno = 0;
const std::string default_name = "";
const int default_continuous = 0;

Sample::Sample(int p, std::string n)
   :phenotype(p), name(n)
{
   if (p != 0 && p != 1)
   {
      continuous_pheno = 1;
   }
   else
   {
      continuous_pheno = 0;
   }
}

Sample::Sample()
   :phenotype(default_pheno), continuous_pheno(default_continuous), name(default_name)
{
}

std::istream& operator>>(std::istream &is, Sample& s)
{
   int phenotype;
   std::string FID, IID;

   is >> FID >> IID >> phenotype;
   if (is)
   {
      s = Sample(phenotype, IID);
   }

   return is;
}
