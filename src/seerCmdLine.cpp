/*
 * File: seerCmdLine.cpp
 *
 * Reads command line input to seer using boost
 * program options
 *
 */

#include "seer.hpp"

namespace po = boost::program_options; // Save some typing

// Use boost::program_options to parse command line input
// This does pretty much all the parameter checking needed
int parseCommandLine (int argc, char *argv[], po::variables_map& vm)
{
   int failed = 0;

   //Required options
   po::options_description required("Required options");
   required.add_options()
    ("kmers,k", po::value<std::string>()->required(), "dsm kmer output file")
    ("pheno,p", po::value<std::string>()->required(), ".pheno metadata");

   // kmds options
   po::options_description covar("Covariate options");
   covar.add_options()
    ("struct", po::value<std::string>(), "mds values from kmds")
    ("covar_file", po::value<std::string>(), "file containing covariates")
    ("covar_list", po::value<std::string>(), "list of columns covariates to use. Format is 1,2q,3 (use q for quantitative)");

   //Optional filtering parameters
   //NB pval cutoffs are strings for display, and are converted to floats later
   po::options_description performance("Performance options");
   performance.add_options()
    ("threads", po::value<int>()->default_value(1), ("number of threads. Suggested: " + std::to_string(std::thread::hardware_concurrency())).c_str());

   //Optional filtering parameters
   //NB pval cutoffs are strings for display, and are converted to floats later
   po::options_description filtering("Filtering options");
   filtering.add_options()
    ("no_filtering", "turn off all filtering and perform tests on all kmers input")
    ("max_length", po::value<long int>()->default_value(max_length_default), "maximum kmer length")
    ("maf", po::value<double>()->default_value(maf_default), "minimum kmer frequency")
    ("min_words", po::value<int>(), "minimum kmer occurences. Overrides --maf")
    ("chisq", po::value<std::string>()->default_value(chisq_default), "p-value threshold for initial chi squared test. Set to 1 to show all")
    ("pval", po::value<std::string>()->default_value(pval_default), "p-value threshold for final logistic test. Set to 1 to show all");

   po::options_description other("Other options");
   other.add_options()
    ("print_samples", "print lists of samples significant kmers were found in")
    ("version", "prints version and exits")
    ("help,h", "full help message");

   po::options_description all;
   all.add(required).add(covar).add(performance).add(filtering).add(other);

   try
   {
      po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

      if (vm.count("help"))
      {
         printHelp(all);
         failed = 1;
      }
      else if (vm.count("version"))
      {
         std::cout << VERSION << std::endl;
         failed = 1;
      }
      else
      {
         po::notify(vm);
         failed = 0;

         // Check input files exist, and can stat
         if (!fileStat(vm["kmers"].as<std::string>()) || !fileStat(vm["pheno"].as<std::string>()))
         {
            failed = 1;
         }
         else if (vm.count("struct") && !fileStat(vm["struct"].as<std::string>()))
         {
            failed = 1;
         }
      }

   }
   catch (po::error& e)
   {
      // Report errors from boost library
      std::cerr << "Error in command line input: " << e.what() << "\n";
      std::cerr << "Run 'seer --help' for full option listing\n\n";
      std::cerr << required << "\n" << other << "\n";

      failed = 1;
   }

   return failed;
}

// Print long help message
void printHelp(po::options_description& help)
{
   std::cerr << help << "\n";
}

