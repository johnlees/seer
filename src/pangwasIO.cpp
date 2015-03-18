/*
 * File: pangwasIO.cpp
 *
 * Controls: reading in pheno to sample object
 * Internal format conversion
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

arma::vec constructVecY(const std::vector<Sample>& samples)
{
   arma::vec y;
   y.zeros(samples.size());

   for (int i = 0; i < samples.size(); ++i)
   {
      y[i] = samples[i].pheno();
   }

   return y;
}

// Vectors x and y for linear relationship y = bx
arma::vec constructVecX(const Kmer& k, const std::vector<Sample>& samples)
{
   arma::vec x;
   x.zeros(samples.size());

   // Names in k and samples are sorted in the same order
   std::vector<std::string>::iterator nameIt = k.occurrence_vector().begin();
   for (int i = 0; i < samples.size(); ++i)
   {
      if(*nameIt == samples[i].iid())
      {
         nameIt++;
         x[i] = 1;
      }
   }

   return x;
}






