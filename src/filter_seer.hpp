/*
 * Header file for filter_seer.cpp
 * Filters significant kmers from seer
 *
 */

// C++ stl includes
#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <iterator>
#include <stdexcept>
#include <math.h>

// Library includes
#include <boost/program_options.hpp>

#include "significant_kmer.hpp"

// Structures
struct cmdOptions
{
   std::string input_file, output_file, sort_field;
   double maf_filter, chi_filter, p_filter, beta_filter;
   bool neg_beta, substr_kmers;
};

// Function prototypes
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
cmdOptions processCmdLine(boost::program_options::variables_map& vm);
double fractionFilter(const std::string& filter_input);
void printHelp(boost::program_options::options_description& help);

