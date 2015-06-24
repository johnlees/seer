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
#include <queue>
#include <future>

// Constants
//    Default options
const int pc_default = 3;
const long int size_default = 1000000;

// Structs
struct distance_element
{
   unsigned int row;
   unsigned int col;
   double distance;
};

// panglossCmdLine headers
int parseCommandLine (int argc, char *argv[], po::variables_map& vm);
void printHelp(po::options_description& help);

// panglossStruct headers
arma::mat metricMDS(const arma::mat& populationMatrix, const int dimensions, const unsigned int threads, const std::string& distances_file = "");
arma::mat dissimiliarityMatrix(const arma::mat& inMat, const unsigned int threads);

distance_element threadDistance(const unsigned int i, const unsigned int j, const arma::rowvec row_1, const arma::rowvec row_2);
double distanceFunction(const arma::rowvec& vec_1, const arma::rowvec& vec_2);

