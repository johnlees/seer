/*
 * File: sample.cpp
 *
 * Helper functions for the sample class
 *
 */

#include "sample.h"

const int default_pheno = 0;
const std::string default_name = "";

Sample::Sample(int p, std::string n)
   :phenotype(p), name(n)
{
   if (p != 0 && p != 1)
   {
      throw std::runtime_error("Sample " + name + " has invalid phenotype " + std::to_string(p) + "\n");
   }
}

Sample::Sample()
   :phenotype(default_pheno), name(default_name)
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
