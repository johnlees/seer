/*
 * File: kmer.cpp
 *
 * Helper functions for the kmer class
 *
 */

#include "kmer.hpp"

// Initialisation
Kmer::Kmer(const std::string& sequence, const std::vector<std::string>& occurrences, const double pvalue, const double lrt_pvalue, const double beta, const double se, const double maf, const double log_likelihood)
   : Significant_kmer(sequence, occurrences, maf, kmer_chi_pvalue_default, pvalue, lrt_pvalue, beta, se, kmer_comment_default), _x_set(0), _log_likelihood(log_likelihood), _use_firth(0)
{
}

// Initialise without calculated information
Kmer::Kmer(const std::string& sequence, const std::vector<std::string>& occurrences)
   : Significant_kmer(sequence, occurrences, kmer_maf_default, kmer_chi_pvalue_default, kmer_pvalue_default, kmer_pvalue_default, kmer_beta_default, kmer_se_default, kmer_comment_default), _x_set(0), _log_likelihood(0), _use_firth(0)
{
}

// Initialise with default info only
Kmer::Kmer()
    : Significant_kmer(kmer_seq_default, kmer_occ_default, kmer_maf_default, kmer_chi_pvalue_default, kmer_pvalue_default, kmer_pvalue_default, kmer_beta_default, kmer_se_default, kmer_comment_default), _x_set(0), _log_likelihood(0), _use_firth(0)
{
}

// Output
std::ostream& operator<<(std::ostream &os, const Kmer& k)
{
   os << std::fixed << std::setprecision(3) << k.sequence() << "\t" << k.maf()
      << "\t" << std::scientific << k.unadj() << "\t" << k.p_val() << "\t" << k.lrt_p_val()
      << "\t" << k.beta() << "\t" << k.se();

   std::vector<double> covariates = k.covar_p();
   for (auto it = covariates.begin(); it != covariates.end(); ++it)
   {
      os << "\t" << *it;
   }

   os << "\t" << k.comments();

   return os;
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
    *
    * OR
    *
    * AAAAAAAAAAAAAAAAAATGCATATTTATCTTAG 5.172314 6925_3#7:9 6823_4#17:26 6871_2#9:8
    *
    */
   if (is)
   {
      std::getline(is, kmer_line);
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
         // Gone one too far in a string w/out | separator
         else
         {
            size_t occurence_pos = s_buffer.find(":");
            if (s_buffer.find(":") != std::string::npos)
            {
               samples.push_back(s_buffer.substr(0, occurence_pos));
               break;
            }
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

   k = Kmer(word, samples);

   return is;
}

// Set the x and maf of the kmer
void Kmer::add_x(const std::unordered_map<std::string,int>& sample_map, const int num_samples)
{
   _x.zeros(num_samples);

   for (auto it = _samples.begin(); it != _samples.end(); ++it)
   {
      auto sample_index_it = sample_map.find(*it);

      if (sample_index_it != sample_map.end())
      {
         _x[sample_index_it->second] = 1;
      }
   }

   // Set object values
   _x_set = 1;
   _maf = (double)num_occurrences()/num_samples;
}

// Add a new comment in
void Kmer::add_comment(const std::string& new_comment)
{
   if (_comment == kmer_comment_default)
   {
      _comment = new_comment;
   }
   else
   {
      _comment += "," + new_comment;
   }
}

int Kmer::num_occurrences() const
{
   int total_occurrences;
   if (_x_set)
   {
      total_occurrences = (int) accu(_x);
   }
   else
   {
      total_occurrences = _samples.size();
   }

   return total_occurrences;
}

