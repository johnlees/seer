/*
 * kmer.hpp
 * Header file for kmer class
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_map>

#define ARMA_DONT_PRINT_ERRORS
#include <armadillo>

#include "significant_kmer.hpp"

const std::string kmer_seq_default = "";
const std::vector<std::string> kmer_occ_default;
const double kmer_pvalue_default = 1;
const double kmer_chi_pvalue_default = 1;
const double kmer_beta_default = 0;
const double kmer_maf_default = 0;
const double kmer_se_default = 0;
const std::string kmer_comment_default = "NA";

class Kmer: public Significant_kmer
{
   public:
      // Initialisation
      Kmer(const std::string& sequence, const std::vector<std::string>& occurrences, const double pvalue, const double beta, const double se, const double maf);
      Kmer(const std::string& sequence, const std::vector<std::string>& occurrences); // Initialise without calculated information
      Kmer(); // defaults

      // nonmodifying operations
      int length() const { return _word.length(); }
      int num_occurrences() const;
      std::string occurrence(int i) const { return _samples[i]; }
      std::vector<std::string> occurrence_vector() const { return _samples; }
      arma::vec get_x() const { return _x; }
      int has_x() const { return _x_set; }

      // Modifying operations
      void add_comment(const std::string& new_comment); // this is defined in kmer.cpp
      void add_x(const std::unordered_map<std::string,int>& sample_map, const int num_samples); // this is defined in kmer.cpp

   private:
      arma::vec _x;
      int _x_set;

};

// Overload output and input operators
std::ostream& operator<<(std::ostream &os, const Kmer& k);

std::istream& operator>>(std::istream &is, Kmer& k);
