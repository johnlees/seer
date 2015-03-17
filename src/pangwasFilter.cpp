/*
 * File: pangwasFilter.cpp
 *
 * Filters kmers that won't be meaningfully associated.
 * Saves time over performing regressions on all, and reduces false positives
 *
 */

#include "pangwas.h"

// Wrapper to all filter functions
int passFilters(const Kmer& k, const std::vector<Sample>& samples)
{
   int passed = 1;

   return passed;
}
