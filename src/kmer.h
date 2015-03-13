/*
 * kmer.h
 * Header file for kmer class
 */

#include <iostream>
#include <string>
#include <vector>

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
      int num_occurences() const { return occurrences.size(); }
      long int get_position() const { return position; }
      double p_val() const { return pvalue; }
      double beta() const { return b; }

      //Modifying operations
      void p_val(double p) { pvalue = p; }
      void beta(double new_b) { b = new_b; }
      long int map(std::string& ref_file);

   private:
      std::string word;
      std::vector<std::string> occurrences;
      double pvalue;
      double b;
      long int position;
};

// Overload output and input operators
std::ostream& operator<<(std::ostream &os, const Kmer& k);

std::istream& operator>>(std::istream &is, Kmer& k);
