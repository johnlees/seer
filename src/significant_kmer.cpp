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

   // Read the line with sequence and p-value, extract only sequence
   std::getline(is, line_in);
   std::stringstream line_stream(line_in);

   line_stream >> sequence;

   // Read the line with samples, extract using the same format as they were
   // written with
   std::getline(is, line_in);
   line_stream.str(line_in);

   std::copy(std::istream_iterator<std::string>(line_stream), std::istream_iterator<std::string>(), std::back_inserter(sample_list));

   // Ensure vector remains sorted on sample name
   std::sort(sample_list.begin(), sample_list.end());
   sk = Significant_kmer(sequence, sample_list);

   return is;
}
