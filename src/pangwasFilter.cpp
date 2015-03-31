/*
 * File: pangwasFilter.cpp
 *
 * Filters kmers that won't be meaningfully associated.
 * Saves time over performing regressions on all, and reduces false positives
 *
 */

#include "pangwas.h"

// Wrapper to all filter functions
int passBasicFilters(const Kmer& k, const int max_length, const int min_words, const int max_words)
{
   int passed = 1;

   //TODO there might be a nicer way to write this, each filter as its own
   //function, then call all functions of type filter on the passed objects

   // Don't test long kmers
   if (k.length() > max_length)
   {
      passed = 0;
   }

   // Impose min words
   // TODO may want to make this more sophisticated, make sure there are at
   // least ten words in each category
   if (passed && k.num_occurrences() < min_words || k.num_occurrences() > max_words)
   {
      passed = 0;
   }

   return passed;
}

int passStatsFilters(const arma::vec& x, const arma::vec& y, double chi_cutoff)
{
   int passed = 1;

   // Contigency table
   //         unaffected affected
   // present a          b
   // absent  c          d
   int a = 0, b = 0, c = 0, d = 0;

   arma::vec::const_iterator j = y.begin();
   for (arma::vec::const_iterator i = x.begin(); i!=x.end(); ++i)
   {
      if (*j == 0) {
         if (*i == 0){
            c++;
         } else {
            a++;
         }
      } else {
         if (*i == 0){
            d++;
         } else {
            b++;
         }
      }
      j++;
   }

   arma::mat::fixed<2, 2> table = {a, b, c, d};
#ifdef PANGWAS_DEBUG
   std::cerr << table << "\n";
#endif

   if (chiTest(table) > chi_cutoff)
   {
      passed = 0;
   }

   return passed;
}
