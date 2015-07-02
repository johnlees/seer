/*
 *
 * pangwas.hpp
 * Header file for pangwas package
 *
 */

// Common headers
//#include "pancommon.hpp"
#include "logitFunction.hpp" // This includes pancommon.hpp

// dlib headers
#include <dlib/optimization.h>

// Constants
//    Default options
const std::string pval_default = "10e-8";
const double convergence_limit = 10e-8;
const unsigned int max_nr_iterations = 1000;

// pangwasCmdLine headers
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);

// pangwasBinaryAssoc headers
void logisticTest(Kmer& k, const arma::vec& y, const unsigned int nr);
void logisticTest(Kmer& k, const arma::vec& y_train, const unsigned int nr, const arma::mat& mds);

void doLogit(Kmer& k, const arma::vec& y_train, const arma::mat& x_train, const unsigned int nr);

regression logisticPval(const arma::vec& y_train, const arma::mat& x_train);
regression newtonRaphson(const arma::vec& y_train, const arma::mat& x_train);
arma::mat varCovarMat(const arma::mat& x, const arma::mat& b);
arma::vec predictLogitProbs(const arma::mat& x, const arma::vec& b);

// pangwasContinuousAssoc headers
void linearTest(Kmer& k, const arma::vec& y_train);
void linearTest(Kmer& k, const arma::vec& y_train, const arma::mat& mds);

void doLinear(Kmer& k, const arma::vec& y_train, const arma::mat& x_train);
arma::vec predictLinearProbs(const arma::mat& x, const arma::vec& b);

