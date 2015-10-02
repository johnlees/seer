/*
 * kmer.h
 * Header file for kmer class
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_map>

#include <armadillo>

const std::string kmer_seq_default = "";
const std::vector<std::string> kmer_occ_default;
const double kmer_pvalue_default = 1;
const double kmer_chi_pvalue_default = 1;
const double kmer_beta_default = 0;
const double kmer_maf_default = 0;

class Kmer
{
   public:
      // Initialisation
      Kmer(std::string sequence, std::vector<std::string> occurrences, double pvalue, double beta, double maf);
      Kmer(std::string sequence, std::vector<std::string> occurrences); // Initialise without calculated information
      Kmer(); // defaults

      // nonmodifying operations
      std::string sequence() const { return word; }
      int length() const { return word.length(); }
      int num_occurrences() const { return occurrences.size(); }
      std::string occurrence(int i) const { return occurrences[i]; }
      std::vector<std::string> occurrence_vector() const { return occurrences; }
      arma::vec get_x() const { return x; }
      int has_x() const { return x_set; }
      double p_val() const { return pvalue; }
      double chi_p_val() const { return chi_pvalue; }
      double beta() const { return b; }
      double get_maf() const { return maf; }

      //Modifying operations
      void p_val(const double _pvalue) { pvalue = _pvalue; }
      void chi_p_val(const double _chi_pvalue) { chi_pvalue = _chi_pvalue; }
      void beta(const double _b) { b = _b; }
      void set_maf(const double _maf) { maf = _maf; }
      void add_x(const std::unordered_map<std::string,int>& sample_map, const int num_samples); // this is defined in kmer.cpp

   private:
      std::string word;
      std::vector<std::string> occurrences;
      arma::vec x;
      int x_set;
      double pvalue;
      double chi_pvalue;
      double b;
      double maf;
};

// Overload output and input operators
std::ostream& operator<<(std::ostream &os, const Kmer& k);

std::istream& operator>>(std::istream &is, Kmer& k);
