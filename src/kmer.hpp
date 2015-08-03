/*
 * kmer.h
 * Header file for kmer class
 */

#include <iostream>
#include <string>
#include <vector>

#include <armadillo>

const std::string kmer_seq_default = "";
const std::vector<std::string> kmer_occ_default;
const double kmer_pvalue_default = 1;
const double kmer_chi_pvalue_default = 1;
const double kmer_beta_default = 0;
const int kmer_position_default = 0;

class Kmer
{
   public:
      // Initialisation
      Kmer(std::string sequence, std::vector<std::string> occurrences, double pvalue, double beta, long int position);
      Kmer(std::string sequence, std::vector<std::string> occurrences); // Initialise without calculated information
      Kmer(); // defaults

      // nonmodifying operations
      std::string sequence() const { return word; }
      int length() const { return word.length(); }
      int num_occurrences() const { return occurrences.size(); }
      std::string occurrence(int i) const { return occurrences[i]; }
      std::vector<std::string> occurrence_vector() const { return occurrences; }
      long int get_position() const { return position; }
      arma::vec get_x() const { return x; }
      int has_x() const { return x_set; }
      double p_val() const { return pvalue; }
      double chi_p_val() const { return chi_pvalue; }
      double beta() const { return b; }

      //Modifying operations
      void p_val(double p) { pvalue = p; }
      void chi_p_val(double _chi_pvalue) { chi_pvalue = _chi_pvalue; }
      void beta(double new_b) { b = new_b; }
      void add_x(arma::vec new_x) { x = new_x; x_set = 1; }
      long int map(std::string& ref_file);

   private:
      std::string word;
      std::vector<std::string> occurrences;
      arma::vec x;
      int x_set;
      double pvalue;
      double chi_pvalue;
      double b;
      long int position;
};

// Overload output and input operators
std::ostream& operator<<(std::ostream &os, const Kmer& k);

std::istream& operator>>(std::istream &is, Kmer& k);
