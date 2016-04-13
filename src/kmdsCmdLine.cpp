/*
 * File: kmdsCmdLine.cpp
 *
 * Reads command line input to kmds
 *
 */

#include "kmds.hpp"

namespace po = boost::program_options; // Save some typing

// Use boost::program_options to parse command line input
// This does pretty much all the parameter checking needed
int parseCommandLine (int argc, char *argv[], po::variables_map& vm)
{
   int failed = 0;

   //Required options
   po::options_description required("Required options");
   required.add_options()
    ("kmers,k", po::value<std::string>(), "dsm kmer output file (not needed if using --mds_concat)")
    ("pheno,p", po::value<std::string>(), ".pheno metadata");

   po::options_description mds("MDS options");
   mds.add_options()
    ("output,o", po::value<std::string>(), "output prefix for new dsm file")
    ("no_mds", "do not perform MDS; output subsampled matrix instead")
    ("write_distances",  "write csv of distance matrix")
    ("mds_concat", po::value<std::string>(), "list of subsampled matrices to use in MDS. Performs only MDS; implies --no_filtering")
    ("pc", po::value<int>()->default_value(pc_default), "number of principal coordinates to output")
    ("size", po::value<long int>()->default_value(size_default), "number of kmers to use in MDS")
    ("threads", po::value<int>()->default_value(1), ("number of threads. Suggested: " + std::to_string(std::thread::hardware_concurrency())).c_str());

   //Optional filtering parameters
   //NB pval cutoffs are strings for display, and are converted to floats later
   po::options_description filtering("Filtering options");
   filtering.add_options()
    ("no_filtering", "turn off all filtering and do not output new kmer file")
    ("max_length", po::value<long int>()->default_value(max_length_default), "maximum kmer length")
    ("maf", po::value<double>()->default_value(maf_default), "minimum kmer frequency")
    ("min_words", po::value<int>(), "minimum kmer occurences. Overrides --maf");

   po::options_description other("Other options");
   other.add_options()
    ("version", "prints version and exits")
    ("help,h", "full help message");

   po::options_description all;
   all.add(required).add(mds).add(filtering).add(other);

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
         if ((vm.count("kmers") && !fileStat(vm["kmers"].as<std::string>())) || (vm.count("pheno") && !fileStat(vm["pheno"].as<std::string>())))
         {
            failed = 1;
         }
      }

   }
   catch (po::error& e)
   {
      // Report errors from boost library
      std::cerr << "Error in command line input: " << e.what() << "\n";
      std::cerr << "Run 'kmds --help' for full option listing\n\n";
      std::cerr << required << "\n" << other << "\n";

      failed = 1;
   }

   return failed;
}

// Print long help message
void printHelp(po::options_description& help)
{
   std::cerr << "kmds" << "\n";
   std::cerr << "\t1) filter and subsample with --no_mds and --size\n";
   std::cerr << "\t2) combine, and do metric multidimensional scaling with --mds_concat\n";

   std::cerr << help << "\n";
}

