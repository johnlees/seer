/*
 * File: pangwasMain.cpp
 *
 * Reads command line input to pangwas, and controls relevant functions
 *
 */

#include "pangwas.h"

namespace po = boost::program_options; // Save some typing

int main (int argc, char *argv[])
{
   // Program description
   std::cerr << "pangwas: pan-genome wide association study using kmers\n";

   po::variables_map vm;

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   if (argc == 1)
   {
      std::cerr << "Usage: pangwas -k dsm.txt -p data.pheno -o output\n\n"
         << "For full option details run pangwas -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, &vm))
   {
      return 1;
   }

   if (vm.count("compression"))
   {
      std::cout << "Compression level was set to " << vm["compression"].as<int>() << ".\n";
   }
   else
   {
      std::cout << "Compression level was not set.\n";
   }
}

// Use boost::program_options to parse command line input
int parseCommandLine (int argc, char *argv[], po::variables_map *vm)
{
   int failed = 0;

   //TODO look at vector string vs. string
   //Required options
   po::options_description required("Required options");
   required.add_options()
    ("help,h", "full help message")
    ("kmers,k", po::value<std::string>()->required(), "dsm kmer output file")
    ("pheno,p", po::value<std::string>()->required(), ".pheno metadata")
    ("output,o", po::value<std::string>()->required(), "output prefix");

   //Optional filtering parameters
   //NB pval cutoffs are strings for display, and are converted to floats later
   po::options_description filtering("Filtering options");
   filtering.add_options()
    ("max_length", po::value<int>(), "maximum kmer length")
    ("maf", po::value<double>()->default_value(maf_default), "minimum kmer frequency")
    ("min_words", po::value<int>(), "minimum kmer occurences. Overrides --maf")
    ("chi2", po::value<std::string>()->default_value(chi2_default), "p-value threshold for initial chi squared test")
    ("pval", po::value<std::string>()->default_value(pval_default), "p-value threshold for final logistic test");

   po::options_description all;
   all.add(required).add(filtering);

   try
   {
      po::store(po::command_line_parser(argc, argv).options(all).run(), *vm);

      if (vm->count("help"))
      {
         printHelp(&all);
         failed = 1;
      }
      else
      {
         po::notify(*vm);
         failed = 0;
      }

   }
   catch (po::error& e)
   {
      std::cerr << "Error in command line input: " << e.what() << "\n";
      std::cerr << "Run 'pangwas --help' for full option listing\n\n";
      std::cerr << required << "\n";

      failed = 1;
   }

   return failed;
}

// Print long help message
void printHelp(po::options_description *help)
{
   std::cerr << *help << "\n";
}
