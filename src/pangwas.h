/*
 *
 * pangwas.h
 * Header file for pangwas package
 *
 */

// C headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <list>
#include <vector>
#include <thread>
#include <exception>
#include <sys/stat.h>

// Boost headers
#include <boost/program_options.hpp>

// Classes
#include "kmer.h"
#include "sample.h"

// Constants
//    Default options
const long int max_length_default = 100;
const double maf_default = 0.01;
const std::string chi2_default = "10e-5";
const std::string pval_default = "10e-8";

// Structs

// Association results return a pvalue and some form of effect size
// Take the slope (beta_1)
struct Assoc
{
   double pvalue;
   double beta;
};

// pangwasMain headers
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);

void printHelp(boost::program_options::options_description& help);
int fileStat(const std::string& filename);

// pangwasIO headers
void readPheno(const std::string& filename, std::vector<Sample>& samples);

// pangwasFilter headers

// pangwasAssoc headers

