/*
 *
 * map_back.hpp
 * Header file for map_back
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <iterator>
#include <vector>
#include <stdexcept>

// Boost headers
#include <boost/program_options.hpp>

// Classes
#include "fasta.hpp"
#include "significant_kmer.hpp"

// Constants

// Structs

// Function headers
// mapMain.cpp
std::vector<Fasta> readSequences(const std::string& reference_file);

// mapCmdLine.cpp
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);

