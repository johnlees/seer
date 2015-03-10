/*
 *
 * pangwas.h
 * Header file for pangwas package
 *
 */

// C headers
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <list>
#include <vector>
#include <thread>
#include <sys/stat.h>

// Boost headers
#include <boost/program_options.hpp>

// Constants
//    Default options
const long int max_length_default = 100;
const double maf_default = 0.01;
const std::string chi2_default = "10e-5";
const std::string pval_default = "10e-8";

// pangwasMain headers
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);

void printHelp(boost::program_options::options_description& help);
int fileStat(const std::string& filename);

