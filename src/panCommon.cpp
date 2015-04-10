/*
 * File: panCommon.cpp
 *
 * Functions common to pangwas and pangloss
 * Program options parsing
 *
 */

#include "pancommon.hpp"

// Parse command line parameters into usable program parameters
cmdOptions verifyCommandLine(boost::program_options::variables_map& vm, const std::vector<Sample>& samples)
{
   cmdOptions verified;

   verified.chi_cutoff = stod(vm["chi2"].as<std::string>());
   verified.max_length = vm["max_length"].as<long int>();

#ifndef NO_THREADS
   if (vm.count("threads") && vm["threads"].as<int>() > 0)
   {
      verified.num_threads = vm["threads"].as<int>();
   }
   else
   {
      verified.num_threads = 1;
   }
#else
   verified.num_threads = 1;
#endif

   if(vm.count("kmers"))
   {
      verified.kmers = vm["kmers"].as<std::string>();
   }

   if(vm.count("output"))
   {
      verified.output = vm["output"].as<std::string>();
   }

   if (vm.count("pval"))
   {
      verified.log_cutoff = stod(vm["pval"].as<std::string>());
   }

   if (vm.count("pc"))
   {
      if (vm["pc"].as<int>() > 0)
      {
         verified.pc = vm["pc"].as<int>();
      }
      else
      {
         badCommand("pc", std::to_string(vm["pc"].as<int>()));
      }
   }

   if (vm.count("size"))
   {
      if (vm["size"].as<long int>() > 0)
      {
         verified.size = vm["size"].as<long int>();
      }
      else
      {
         badCommand("size", std::to_string(vm["size"].as<long int>()));
      }
   }

   verified.filter = 1;
   if (vm.count("no_filtering"))
   {
      verified.filter = 0;
   }

   // Error check cmd line options
   verified.min_words = 0;
   if (vm.count("min_words"))
   {
      int min_words_in = vm["min_words"].as<int>();
      if (min_words_in >= 0)
      {
         verified.min_words = min_words_in;
      }
      else
      {
         badCommand("min_words", std::to_string(min_words_in));
      }
   }
   else
   {
      double maf_in = vm["maf"].as<double>();
      if (maf_in >= 0)
      {
         verified.min_words = static_cast<unsigned int>(samples.size() * maf_in);
      }
      else
      {
         badCommand("maf", std::to_string(maf_in));
      }
   }

   if (verified.min_words > samples.size())
   {
      badCommand("min_words/maf", std::to_string(verified.min_words));
   }
   verified.max_words = samples.size() - verified.min_words;

   return verified;
}
