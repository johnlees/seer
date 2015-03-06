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

// Boost headers
#include <boost/program_options.hpp>

// pangwasMain headers
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map *vm);

