/*
 * filterCmdLine.cpp
 * Parses cmd line options for filter_seer
 *
 */

#include "filter_seer.hpp"

namespace po = boost::program_options; // Save some typing

// Parse command line options using boost program options
int parseCommandLine (int argc, char *argv[], po::variables_map& vm)
{
   int failed = 0;

   //Required options
   po::options_description required("Required options");
   required.add_options()
    ("kmers,k", po::value<std::string>()->required(), "file of output from seer");

   po::options_description filters("Filtering options");
   filters.add_options()
    ("chisq", po::value<std::string>(), "minimum unadjusted p-value to output")
    ("pval", po::value<std::string>(), "minimum adjusted p-value to output")
    ("maf", po::value<std::string>(), "minimum maf/max 1-maf to output")
    ("beta", po::value<std::string>(), "minimum |beta| to output")
    ("substr", "remove smaller kmers completely represented elsewhere")
    ("pos_beta", "output positive effect sizes only");

   po::options_description other("Other options");
   other.add_options()
    ("sort,s", po::value<std::string>(), "field to sort on: chisq, pval, maf or beta")
    ("help,h", "full help message");

   po::options_description all;
   all.add(required).add(filters).add(other);

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

// Makes sure command line options are valid. Casts into correct variables
cmdOptions processCmdLine(po::variables_map& vm)
{
   cmdOptions processed_options;

   processed_options.input_file = vm["kmers"].as<std::string>();

   if (vm.count("chisq"))
   {
      processed_options.chi_filter = fractionFilter(vm["chisq"].as<std::string>());
   }
   else
   {
      processed_options.chi_filter = 0;
   }

   if (vm.count("pval"))
   {
      processed_options.p_filter = fractionFilter(vm["pval"].as<std::string>());
   }
   else
   {
      processed_options.p_filter = 0;
   }

   if (vm.count("maf"))
   {
      processed_options.maf_filter = fractionFilter(vm["maf"].as<std::string>());
      if (processed_options.maf_filter > 0.5)
      {
         processed_options.maf_filter = 1 - processed_options.maf_filter;
      }
   }
   else
   {
      processed_options.maf_filter = 0;
   }

   if (vm.count("beta"))
   {
      processed_options.beta_filter = stod(vm["beta"].as<std::string>());
   }
   else
   {
      processed_options.beta_filter = 0;
   }

   if (vm.count("pos_beta"))
   {
      processed_options.neg_beta = 0;
   }
   else
   {
      processed_options.neg_beta = 1;
   }

   if (vm.count("substr"))
   {
      processed_options.substr_kmers = 0;
   }
   else
   {
      processed_options.substr_kmers = 1;
   }

   if (vm.count("sort"))
   {
      processed_options.sort_field = vm["sort"].as<std::string>();
   }
   else
   {
      processed_options.sort_field = "";
   }

   return processed_options;
}

// Turns a string into a double between 0 and 1
double fractionFilter(const std::string& filter_input)
{
   double parsed_filter = stod(filter_input);
   if (parsed_filter < 0 || parsed_filter > 1)
   {
      parsed_filter = 0;
   }

   return parsed_filter;
}

// Print long help message
void printHelp(po::options_description& help)
{
   std::cerr << help << "\n";
}
