/*
 *
 * pangwas.h
 * Header file for pangwas package
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

// Boost headers
#include <boost/program_options.hpp>
#include <boost/math/distributions/normal.hpp>

// gzstream headers
#include <gzstream.h>

// Armadillo/mlpack headers
#include <mlpack/core.hpp> // this includes armadillo
#include <mlpack/methods/logistic_regression/logistic_regression.hpp>

// Classes
#include "kmer.hpp"
#include "sample.hpp"

// Constants
//    Default options
const long int max_length_default = 100;
const double maf_default = 0.01;
const std::string chi2_default = "10e-5";
const std::string pval_default = "10e-8";

// pangwasCmdLine headers
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);

// pangwasIO headers
void readPheno(const std::string& filename, std::vector<Sample>& samples);
void openDsmFile(igzstream& dsm_file, const std::string file_name);

arma::vec constructVecY(const std::vector<Sample>& samples);
arma::vec constructVecX(const Kmer& k, const std::vector<Sample>& samples);

int fileStat(const std::string& filename);

// pangwasFilter headers
int passBasicFilters(const Kmer& k, const int max_length, const int min_words, const int max_words);
int passStatsFilters(const arma::vec& x, const arma::vec& y, double chi_cutoff);

// pangwasAssoc headers
void logisticTest(Kmer& k, const arma::vec& y);
arma::mat varCovarMat(const arma::mat& x, const arma::mat& b);
arma::vec predictLogitProbs(const arma::mat& x, const arma::vec& b);

double chiTest(arma::mat& table);
double normalPval(double testStatistic);

//pangwasErr headers
void badCommand(const std::string& command, const std::string& value);
