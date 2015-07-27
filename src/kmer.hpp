/*
 * kmer.h
 * Header file for kmer class
 */

#include <iostream>
#include <string>
#include <vector>

#include <armadillo>

#include <boost/multiprecision/mpfr.hpp>

class Kmer
{
   public:
      // Initialisation
      Kmer(std::string sequence, std::vector<std::string> occurrences, boost::multiprecision::mpfr_float_500 pvalue, boost::multiprecision::mpfr_float_500 chi_pvalue, double beta, long int position);
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
      boost::multiprecision::mpfr_float_500 p_val() const { return pvalue; }
      boost::multiprecision::mpfr_float_500 chi_p_val() const { return chi_pvalue; }
      double beta() const { return b; }

      //Modifying operations
      void p_val(boost::multiprecision::mpfr_float_500 p) { pvalue = p; }
      void chi_p_val(boost::multiprecision::mpfr_float_500 _chi_pvalue) { chi_pvalue = _chi_pvalue; }
      void beta(double new_b) { b = new_b; }
      void add_x(arma::vec new_x) { x = new_x; x_set = 1; }
      long int map(std::string& ref_file);

   private:
      std::string word;
      std::vector<std::string> occurrences;
      arma::vec x;
      int x_set;
      boost::multiprecision::mpfr_float_500 chi_pvalue;
      boost::multiprecision::mpfr_float_500 pvalue;
      double b;
      long int position;
};

// Overload output and input operators
std::ostream& operator<<(std::ostream &os, const Kmer& k);

std::istream& operator>>(std::istream &is, Kmer& k);
