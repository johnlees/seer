/*
 *
 * seercommon.hpp
 * Header file for seercommon
 * Shared functions between seer and kmds
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <chrono>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <thread>
#include <exception>
#include <sys/stat.h>
#include <regex>

// gzstream headers
#include <gzstream.h>

// Boost headers
#include <boost/program_options.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>

// Armadillo/dlib headers
#include <armadillo>
#include <dlib/matrix.h>

// Classes
#include "kmer.hpp"
#include "sample.hpp"

// Constants
//    Default options
const double maf_default = 0.01;
const long int max_length_default = 100;
const std::string chi2_default = "10e-5";

typedef dlib::matrix<double,0,1> column_vector;

// Structs
struct cmdOptions
{
   double log_cutoff;
   double chi_cutoff;

   long int max_length;
   long int size;
   int filter;
   int positive;
   int pc;
   int print_samples;
   int write_distances;
   unsigned int num_threads;
   size_t min_words;
   size_t max_words;

   std::string pheno;
   std::string kmers;
   std::string output;
};

struct regression
{
   double p_val;
   double beta;
};

// Function headers
//    seerCommon.cpp
cmdOptions verifyCommandLine(boost::program_options::variables_map& vm, const std::vector<Sample>& samples);
void verifyMDSOptions(cmdOptions& verified, boost::program_options::variables_map& vm);

arma::vec dlib_to_arma(const column_vector& dlib_vec);
column_vector arma_to_dlib(const arma::vec& arma_vec);

int continuousPhenotype (const std::vector<Sample>& sample_list);

// seerErr headers
void badCommand(const std::string& command, const std::string& value);

// seerIO headers
void readPheno(const std::string& filename, std::vector<Sample>& samples, std::unordered_map<std::string,int>& sample_map);
void openDsmFile(igzstream& dsm_file, const std::string& file_name);

arma::vec constructVecY(const std::vector<Sample>& samples);
arma::vec constructVecX(const Kmer& k, const std::vector<Sample>& samples);

void writeMDS(const std::string& file_name, const arma::mat& MDS);
void writeDistances(const std::string& file_name, const arma::mat& distances);
arma::mat readMDS(const std::string& file_name);
arma::mat readMDSList(const std::string& filename);

int fileStat(const std::string& filename);

// seerFilter headers
int passFilters(const cmdOptions& filterOptions, Kmer& k, const std::vector<Sample>& samples, const std::unordered_map<std::string,int>& sample_map, const arma::vec& y, const int continuous_phenotype);
int passBasicFilters(const Kmer& k, const int max_length, const int min_words, const int max_words);
int passStatsFilters(const arma::vec& x, const arma::vec& y, const double chi_cutoff, const int continuous_phenotype, const int positive_only);

// seerChiFilter headers
double chiTest(const arma::vec& x, const arma::vec& y);
double welchTwoSamplet(const arma::vec& x, const arma::vec& y);
double normalPval(double testStatistic);

