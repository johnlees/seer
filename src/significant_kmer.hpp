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

// kmers have a sequence, and a list of samples they appear in
class Significant_kmer
{
   public:
      // Initialisation
      Significant_kmer();
      Significant_kmer(const std::string& word, const std::vector<std::string>& samples, const double maf, const double unadj_p, const double adj_p, const double beta);

      // nonmodifying operations
      std::vector<std::string> samples_found() const { return _samples; }
      std::string sequence() const { return _word; }
      double maf() const { return _maf; }
      double unadj() const { return _unadj_p; }
      double p_val() const { return _adj_p; }
      double beta() const { return _beta; }

   private:
      std::string _word;
      std::vector<std::string> _samples;

      double _maf;
      double _unadj_p;
      double _adj_p;
      double _beta;
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

