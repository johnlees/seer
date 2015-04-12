/*
 *
 * pangwas.hpp
 * Header file for pangwas package
 *
 */

// Common headers
#include "pancommon.hpp"

// Constants
//    Default options
const std::string pval_default = "10e-8";

// pangwasCmdLine headers
int parseCommandLine (int argc, char *argv[], boost::program_options::variables_map& vm);
void printHelp(boost::program_options::options_description& help);

// pangwasAssoc headers
void logisticTest(Kmer& k, const arma::vec& y);
void logisticTest(Kmer& k, const arma::vec& y_train, const arma::mat& mds);

regression logisticPval(const arma::vec& y_train, const arma::mat& x_train);
arma::mat varCovarMat(const arma::mat& x, const arma::mat& b);
arma::vec predictLogitProbs(const arma::mat& x, const arma::vec& b);

