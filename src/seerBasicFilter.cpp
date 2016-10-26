/*
 * File: seerFilter.cpp
 *
 * Can be filtered in kmds or seer
 * Filters kmers that won't be meaningfully associated.
 * Saves time over performing regressions on all, and reduces false positives
 *
 */

#include "seercommon.hpp"

// k-mer length and frequency
int passBasicFilters(const cmdOptions& filterOptions, const Kmer& k)
{
   int passed = 1;

   // Don't test long kmers
   if (k.length() > filterOptions.max_length)
   {
      passed = 0;
   }

   // Impose min words
   // TODO may want to make this more sophisticated, make sure there are at
   // least ten words in each category
   if (passed && (k.num_occurrences() < filterOptions.min_words || k.num_occurrences() > filterOptions.max_words))
   {
      passed = 0;
   }

   return passed;
}

