/*
 * File: significant_kmer.cpp
 *
 * Helper functions for the significant kmer class
 * Reads input into class
 *
 */

#include "significant_kmer.hpp"

Significant_kmer::Significant_kmer()
{
}

Significant_kmer::Significant_kmer(const std::string& _word, const std::vector<std::string>& _samples)
   :word(_word), samples(_samples)
{
}

// Fills significant kmer object from pangwas output file
// Sample vector is returned sorted
std::istream& operator>>(std::istream &is, Significant_kmer& sk)
{
   std::string sequence, sample, line_in = "";
   std::vector<std::string> sample_list;

   std::getline(is, line_in);
   std::stringstream line_stream(line_in);

   line_stream >> sequence;

   std::getline(is, line_in);
   line_stream.str(line_in);

   while (!line_stream.eof())
   {
      line_stream >> sample;
      sample_list.push_back(sample);
   }

   std::sort(sample_list.begin(), sample_list.end());
   sk = Significant_kmer(sequence, sample_list);

   return is;
}
