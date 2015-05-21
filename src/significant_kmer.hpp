/*
 * significant_kmer.hpp
 * Header file for significant_kmer class
 *
 */

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

// kmers have a sequence, and a list of samples they appear in
class Significant_kmer
{
   public:
      // Initialisation
      Significant_kmer();
      Significant_kmer(const std::string& _word, const std::vector<std::string>& _samples);

      // nonmodifying operations
      std::vector<std::string> samples_found() const { return samples; }
      std::string sequence() const { return word; }

   private:
      std::string word;
      std::vector<std::string> samples;
};

// Overload input operator
std::istream& operator>>(std::istream &is, Significant_kmer& sk);
