/*
 *
 * map_back.hpp
 * Header file for map_back
 *
 */

// C/C++/C++11 headers
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iterator>
#include <vector>
#include <functional>
#include <stdexcept>
#include <future>
#include <thread>
#include <chrono>
#include <list>
#include <assert.h>

// Boost headers
#include <boost/program_options.hpp>

// Classes
#include "fasta.hpp"
#include "significant_kmer.hpp"

// Constants
const std::string VERSION = "1.1.1";
const int thread_wait = 100; //time to wait for each thread in ms

// Structs

// Function headers
// mapMain.cpp
std::vector<Fasta> readSequences(const std::string& reference_file);
void waitForThreads(std::list<std::future<void>>& thread_list, const size_t leave_running);

// mapCmdLine.cpp
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);

