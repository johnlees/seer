/*
 * File: kmer.cpp
 *
 * Helper functions for the kmer class
 *
 */

#include "kmer.h"

const std::string seq_default = "";
const std::vector<std::string> occ_default;
const double pvalue_default = 1;
const double beta_default = 0;
const int position_default = 0;

// Initialisation
Kmer::Kmer(std::string sequence, std::vector<std::string> occurrences, double pvalue, double beta, long int position)
   :word(sequence), occurrences(occurrences), pvalue(pvalue), b(beta), position(position)
{
}

// Initialise without calculated information
Kmer::Kmer(std::string sequence, std::vector<std::string> occurrences)
   :word(sequence), occurrences(occurrences), pvalue(pvalue_default), b(beta_default), position(position_default)
{
}

// Initialise with default info only
Kmer::Kmer()
   :word(seq_default), occurrences(occ_default), pvalue(pvalue_default), b(beta_default), position(position_default)
{
}

// Map kmer to a reference file
long int Kmer::map(std::string& ref_file)
{
   throw std::logic_error("kmer position mapping not yet implemented");
}

// Output
std::ostream& operator<<(std::ostream &os, const Kmer& k)
{
   return os << k.sequence() << "\t" << std::scientific << std::to_string(k.p_val()) << "\t" << k.beta() << "\n";
}

// Input
std::istream& operator>>(std::istream &is, Kmer& k)
{
   std::vector<std::string> samples;
   std::string word, s_buffer, kmer_line;
   char c_buffer;

   /*
    * Example dsm file line AAAAAAAAAAAAAAAAAATGCATATTTATCTTAG 5.172314 0.175087 100 0 100
    * 0.164875 100 | 6925_3#7:9 6823_4#17:26 6871_2#9:8
    */
   if (is)
   {
      getline(is, kmer_line);
      std::stringstream line(kmer_line);

      // First field is the kmer
      line >> word;

      // Skip entropy fields
      while (line >> s_buffer)
      {
         if (s_buffer == "|")
         {
            break;
         }
      }
      s_buffer = "";

      while (line >> c_buffer)
      {
         // Found end of sample. Read forward to next one
         if (c_buffer == ':')
         {
            samples.push_back(s_buffer);
            line >> s_buffer;

            s_buffer = "";
         }
         else
         {
            s_buffer += c_buffer;
         }
      }
   }

   // Keep samples sorted, for easy conversion into a vector
   std::sort(samples.begin(), samples.end());

   k = Kmer(word, samples);

   return is;
}
