/*
 *
 * pancommon.hpp
 * Header file for pancommon
 * Shared functions between pangwas and pangloss
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <list>
#include <vector>
#include <thread>
#include <exception>
#include <sys/stat.h>
#include <regex>

// gzstream headers
#include <gzstream.h>

// Boost headers
#include <boost/program_options.hpp>
#include <boost/math/distributions/normal.hpp>

// Armadillo/mlpack headers
#include <mlpack/core.hpp> // this includes armadillo
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>

// Classes
#include "kmer.hpp"
#include "sample.hpp"

// Constants
//    Default options
const double maf_default = 0.01;
const long int max_length_default = 100;
const std::string chi2_default = "10e-5";

// Structs
struct cmdOptions
{
   double log_cutoff, chi_cutoff;
   long int max_length, size;
   int filter, pc;
   size_t min_words, max_words;
   std::string pheno, kmers, output;
   unsigned int num_threads;
};

struct regression
{
   double p_val;
   double beta;
};

// Function headers
//    panCommon.cpp
cmdOptions verifyCommandLine(boost::program_options::variables_map& vm, const std::vector<Sample>& samples);

// panErr headers
void badCommand(const std::string& command, const std::string& value);

// panIO headers
void readPheno(const std::string& filename, std::vector<Sample>& samples);
void openDsmFile(igzstream& dsm_file, const std::string& file_name);

arma::vec constructVecY(const std::vector<Sample>& samples);
arma::vec constructVecX(const Kmer& k, const std::vector<Sample>& samples);

void writeMDS(const std::string& file_name, const arma::mat& MDS);
arma::mat readMDS(const std::string& file_name);

int fileStat(const std::string& filename);

// panFilter headers
int passFilters(const cmdOptions& filterOptions, Kmer& k, const std::vector<Sample>& samples, const arma::vec& y);
int passBasicFilters(const Kmer& k, const int max_length, const int min_words, const int max_words);
int passStatsFilters(const arma::vec& x, const arma::vec& y, double chi_cutoff);

// panChiFilter headers
double chiTest(arma::mat& table);
double normalPval(double testStatistic);

