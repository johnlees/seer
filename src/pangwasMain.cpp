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

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   po::variables_map vm;
   if (argc == 1)
   {
      std::cerr << "Usage: pangwas -k dsm.txt -p data.pheno -o output\n\n"
         << "For full option details run pangwas -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, vm))
   {
      return 1;
   }

   // Open .pheno file, parse into vector of samples
   std::vector<Sample> samples;
   readPheno(vm["pheno"].as<std::string>(), samples);

   // Open the dsm kmer ifstream
   std::ifstream ist(vm["kmers"].as<std::string>().c_str());
   if (!ist)
   {
      throw std::runtime_error("Could not open kmer file " + vm["kmers"].as<std::string>() + "\n");
   }

   // Parse a set of dsm lines
   std::vector<Kmer> kmer_lines(vm["threads"].as<int>());
   Kmer k;

   for (int i = 0; i<vm["threads"].as<int>(); ++i)
   {
      if (!ist)
      {
         ist >> k;
         kmer_lines.push_back(k);
      }
   }

   // Thread from here...
   // Filter

   // Test

   // ...to here
   // Print in order when all threads complete

   }

// Use boost::program_options to parse command line input
// This does pretty much all the parameter checking needed
int parseCommandLine (int argc, char *argv[], po::variables_map& vm)
{
   int failed = 0;

   //Required options
   po::options_description required("Required options");
   required.add_options()
    ("kmers,k", po::value<std::string>()->required(), "dsm kmer output file")
    ("pheno,p", po::value<std::string>()->required(), ".pheno metadata")
    ("output,o", po::value<std::string>()->required(), "output prefix");

   po::options_description other("Other options");
   other.add_options()
    ("help,h", "full help message");

   //Optional filtering parameters
   //NB pval cutoffs are strings for display, and are converted to floats later
   po::options_description performance("Performance options");
   performance.add_options()
    ("threads", po::value<int>()->default_value(1), "number of threads");

   //Optional filtering parameters
   //NB pval cutoffs are strings for display, and are converted to floats later
   po::options_description filtering("Filtering options");
   filtering.add_options()
    ("max_length", po::value<int>(), "maximum kmer length")
    ("maf", po::value<double>()->default_value(maf_default), "minimum kmer frequency")
    ("min_words", po::value<int>(), "minimum kmer occurences. Overrides --maf")
    ("chi2", po::value<std::string>()->default_value(chi2_default), "p-value threshold for initial chi squared test. Set to zero to show all")
    ("pval", po::value<std::string>()->default_value(pval_default), "p-value threshold for final logistic test. Set to zero to show all");

   po::options_description all;
   all.add(required).add(performance).add(filtering).add(other);

   try
   {
      po::store(po::command_line_parser(argc, argv).options(all).run(), vm);

      if (vm.count("help"))
      {
         printHelp(all);
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
      }

   }
   catch (po::error& e)
   {
      // Report errors from boost library
      std::cerr << "Error in command line input: " << e.what() << "\n";
      std::cerr << "Run 'pangwas --help' for full option listing\n\n";
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

// Check for file existence
int fileStat(const std::string& filename)
{
   struct stat buffer;
   int success = 1;

   if (stat (filename.c_str(), &buffer) != 0)
   {
      std::cerr << "Can't stat input file: " << filename << "\n";

      success = 0;
   }

   return success;
}
