/*
 *
 * kmds.hpp
 * Header file for kmds
 *
 */

// Common headers
#include "seercommon.hpp"

// C/C++ std headers
#include <random>
#include <queue>
#include <future>

// Constants
//    Default options
const int pc_default = 10;
const long int size_default = 1000000;

// Structs
struct DistanceElement
{
   unsigned int row;
   unsigned int col;
   double distance;
};

// kmdsCmdLine headers
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);

// kmdsStruct headers
arma::mat metricMDS(const arma::mat& populationMatrix, const int dimensions, const unsigned int threads, const std::string& distances_file = "");
arma::mat dissimiliarityMatrix(const arma::mat& inMat, const unsigned int threads);

std::vector<DistanceElement> threadDistance(std::vector<DistanceElement> element_list, const arma::mat& rectangular_matrix);
double distanceFunction(const arma::rowvec& vec_1, const arma::rowvec& vec_2);

