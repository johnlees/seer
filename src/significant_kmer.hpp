/*
 * significant_kmer.hpp
 * Header file for significant_kmer class
 *
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <limits>

const int default_covars = 0;
const int standard_cols = 8;

// kmers have a sequence, and a list of samples they appear in
class Significant_kmer
{
   public:
      // Initialisation
      Significant_kmer();
      Significant_kmer(const int num_covars);
      Significant_kmer(const std::string& word, const std::vector<std::string>& samples, const double maf, const double unadj_p, const double adj_p, const double lrt_p, const double beta, const double se, const std::string& comments);

      // nonmodifying operations
      long int line_number() const { return _line_nr; }
      std::vector<std::string> samples_found() const { return _samples; }
      std::string sequence() const { return _word; }
      double maf() const { return _maf; }
      double unadj() const { return _unadj_p; }
      double p_val() const { return _adj_p; }
      double lrt_p_val() const { return _adj_lrt_p; }
      double beta() const { return _beta; }
      double se() const { return _se; }
      std::string comments() const { return _comment; }
      std::vector<double> covar_p() const { return _covar_p; }
      unsigned int num_covars() const; // Defined in significant_kmer.cpp
      std::string rev_comp() const; // Defined in significant_kmer.cpp

      // Modifying operations
      void set_line_nr(const long int line_nr) { _line_nr = line_nr; }
      void p_val(const double pvalue) { _adj_p = pvalue; }
      void lrt_p_val(const double lrt_pvalue) { _adj_lrt_p = lrt_pvalue; }
      void unadj_p_val(const double chi_pvalue) { _unadj_p = chi_pvalue; }
      void beta(const double b) { _beta = b; }
      void standard_error(const double se) { _se = se; }
      void set_maf(const double maf) { _maf = maf; }
      void add_covar_p(const double new_covar_p) { _covar_p.push_back(new_covar_p); }


   protected:
      long int _line_nr;

      std::string _word;
      std::vector<std::string> _samples;

      double _maf;
      double _unadj_p;
      double _adj_p;
      double _adj_lrt_p;
      double _beta;
      double _se;
      std::vector<double> _covar_p;
      std::string _comment;

   private:
      int _num_covars;
};

class sortSigKmer
{
   public:
      sortSigKmer();
      sortSigKmer(const std::string& sort_field);

      void addSortField(const std::string& sort_field);

      // Overload for sort
      bool operator() (const Significant_kmer& sk1, const Significant_kmer& sk2) const;
      explicit operator bool() const { return (_sort_field > 0); };

   private:
      int _sort_field;
};

// Overload input and output operators
std::istream& operator>>(std::istream &is, Significant_kmer& sk);

std::ostream& operator<<(std::ostream &os, const Significant_kmer& k);

// Function for reading header to find number of covariate fields
int parseHeader(const std::string& header_line);

