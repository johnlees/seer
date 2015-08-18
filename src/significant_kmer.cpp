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

   // Read the line convert to stringstream, extract only sequence
   std::getline(is, line_in);
   std::stringstream line_stream(line_in);

   line_stream >> sequence;

   // Ignore fields with p-values, betas etc. Number to ignore set in header
   // file. Add one to this as first tab will be matched by ignore
   for (unsigned int i = 0; i < ignored_fields + 1; ++i)
   {
      line_stream.ignore(std::numeric_limits<std::streamsize>::max(), '\t');
   }

   // Read remainder of line, which is sample names, in the same way they were
   // written
   std::copy(std::istream_iterator<std::string>(line_stream), std::istream_iterator<std::string>(), std::back_inserter(sample_list));

   // Ensure vector remains sorted on sample name
   std::sort(sample_list.begin(), sample_list.end());
   sk = Significant_kmer(sequence, sample_list);

   return is;
}
