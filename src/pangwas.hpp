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

