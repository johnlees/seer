/*
 *
 * pangloss.hpp
 * Header file for pangloss package
 *
 */

// Common headers
#include "pancommon.hpp"

// C/C++ std headers
#include <random>

// Constants
//    Default options
const int pc_default = 3;
const long int size_default = 1000000;

// panglossCmdLine headers
int parseCommandLine (int argc, char *argv[], po::variables_map& vm);
void printHelp(po::options_description& help);

// panglossStruct headers
arma::mat metricMDS(const arma::mat& populationMatrix, const int dimensions);
arma::mat dissimiliarityMatrix(const arma::mat& inMat);

