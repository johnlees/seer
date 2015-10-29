/*
 *
 * seer.hpp
 * Header file for seer package
 *
 */

// Common headers
//#include "seercommon.hpp"
#include "logitFunction.hpp" // This includes seercommon.hpp

// dlib headers
#include <dlib/optimization.h>

// Constants
//    Default options
const std::string pval_default = "10e-8";
const double convergence_limit = 10e-8;
const unsigned int max_nr_iterations = 1000;

// seerCmdLine headers
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);

// seerBinaryAssoc headers
void logisticTest(Kmer& k, const arma::vec& y);
void logisticTest(Kmer& k, const arma::vec& y_train, const arma::mat& mds);

void doLogit(Kmer& k, const arma::vec& y_train, const arma::mat& x_train);
void newtonRaphson(Kmer& k, const arma::vec& y_train, const arma::mat& x_design, const bool firth = 0);

arma::mat varCovarMat(const arma::mat& x, const arma::mat& b);
arma::vec predictLogitProbs(const arma::mat& x, const arma::vec& b);

// seerContinuousAssoc headers
void linearTest(Kmer& k, const arma::vec& y_train);
void linearTest(Kmer& k, const arma::vec& y_train, const arma::mat& mds);

void doLinear(Kmer& k, const arma::vec& y_train, const arma::mat& x_train);
arma::vec predictLinearProbs(const arma::mat& x, const arma::vec& b);

