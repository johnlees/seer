/*
 * combineInit.cpp
 * Functions to initialise inpit for combineKmers
 *
 */

#include "combineKmers.hpp"

std::vector<std::tuple<std::string, std::string> > readSamples(const std::string& sample_file)
{
   // Check input file can be opened, and so so
   std::ifstream ist(sample_file.c_str());

   if (!ist)
   {
      throw std::runtime_error("Could not open sample file " + sample_file + "\n");
   }

   // Structure to store in
   std::vector<std::tuple<std::string, std::string> > samples;
   while(ist)
   {
      std::string sample, file;
      ist >> sample >> file;

      // Don't add after final line
      if (ist)
      {
         samples.push_back(std::make_tuple(sample, file));
      }
   }

   return samples;
}

size_t checkMin(const size_t num_samples, const int input_min_samples)
{
   size_t return_min = input_min_samples;

   if (input_min_samples < 1 || return_min > num_samples)
   {
      return_min = 1;
   }

   return return_min;
}

