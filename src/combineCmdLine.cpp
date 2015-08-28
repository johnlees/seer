/*
 * combineCmdLine.cpp
 * Parses cmd line options for combineKmers
 *
 */

#include "combineKmers.hpp"

namespace po = boost::program_options; // Save some typing

// Parse command line options using boost program options
int parseCommandLine (int argc, char *argv[], po::variables_map& vm)
{
   int failed = 0;

   //Required options
   po::options_description required("Required options");
   required.add_options()
    ("samples,r", po::value<std::string>()->required(), "file with tab separated sample name and kmer file")
    ("output,o", po::value<std::string>()->required(), "output file prefix");

   po::options_description other("Other options");
   other.add_options()
    ("min_samples", po::value<int>()->default_value(1), "minimum number of samples kmer must occur in to be printed")
    ("help,h", "full help message");

   po::options_description all;
   all.add(required).add(other);

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
      }

   }
   catch (po::error& e)
   {
      // Report errors from boost library
      std::cerr << "Error in command line input: " << e.what() << "\n";
      std::cerr << "Run 'combineKmers --help' for full option listing\n\n";
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
