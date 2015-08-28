/*
 * Header file for combineKmers.cpp
 * Takes union of dsk counted kmers
 *
 */

// C++ stl includes
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <iterator>

// Library includes
#include <boost/program_options.hpp>
#include <gzstream.h>

// Function prototypes
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);

std::vector<std::tuple<std::string, std::string> > readSamples(const std::string& sample_file);
size_t checkMin(const size_t num_samples, const int input_min_samples);

