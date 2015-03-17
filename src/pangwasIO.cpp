/*
 * File: pangwasIO.cpp
 *
 * Controls: reading in pheno to sample object
 *
 * NB: some input and output used overridden operators in
 * <class>.cpp files
 *
 */
#include "pangwas.h"

/* .pheno files
 * space separated:
 * FID, IID, phenotype (1 control, 2 case)
 */
void readPheno(const std::string& filename, std::vector<Sample>& samples)
{
   std::ifstream ist(filename.c_str());

   if (!ist)
   {
      throw std::runtime_error("Could not open pheno file " + filename + "\n");
   }

   Sample s;
   while (ist >> s)
   {
      samples.push_back(s);
   }

   // Keep sorted in the same order as the kmers
   std::sort(samples.begin(), samples.end(), Sample::compareSamples);
}






