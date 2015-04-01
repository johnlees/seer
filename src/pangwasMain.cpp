/*
 * File: pangwasMain.cpp
 *
 * Reads command line input to pangwas, and controls relevant functions
 *
 */

#include "pangwas.h"

int main (int argc, char *argv[])
{
   // Program description
   std::cerr << "pangwas: pan-genome wide association study using kmers\n";

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   boost::program_options::variables_map vm;
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

   // Read cmd line options
   double log_cutoff = stod(vm["pval"].as<std::string>());
   double chi_cutoff = stod(vm["chi2"].as<std::string>());
   int max_length = vm["max_length"].as<long int>();
   unsigned int num_threads = 1;

#ifndef NO_THREADS
   if (vm["threads"].as<int>() >= 1)
   {
      num_threads = vm["threads"].as<int>();
   }
#endif

   // Open .pheno file, parse into vector of samples
   std::vector<Sample> samples;
   readPheno(vm["pheno"].as<std::string>(), samples);
   arma::vec y = constructVecY(samples);

   // Error check cmd line options
   size_t min_words = 0;
   if (vm.count("min_words"))
   {
      min_words = vm["min_words"].as<int>();
   }
   else
   {
      min_words = static_cast<int>(samples.size() * vm["maf"].as<double>());
   }

   if (min_words < 0 || min_words > samples.size())
   {
      throw std::runtime_error("Bad min_words/maf specified: " + std::to_string(min_words) + "\n");
   }
   size_t max_words = samples.size() - min_words;

   // Open the dsm kmer ifstream, and read through the whole thing
   std::ifstream kmer_file(vm["kmers"].as<std::string>().c_str());
   if (!kmer_file)
   {
      throw std::runtime_error("Could not open kmer file " + vm["kmers"].as<std::string>() + "\n");
   }

   while (kmer_file)
   {
      // Parse a set of dsm lines
      // TODO this could also be threaded?
      std::vector<Kmer> kmer_lines;
      kmer_lines.reserve(num_threads);

      Kmer k;
      while(kmer_lines.size() < num_threads && kmer_file)
      {
         if (kmer_file)
         {
            kmer_file >> k;

            // apply filters here
            if (passBasicFilters(k, max_length, min_words, max_words))
            {
               // Don't bother with this if not running stats tests
               int pass = 1;
               arma::vec x = constructVecX(k, samples);

               try  // Some chi^2 tests may diverge - proceed anyway for now
               {
                  pass = passStatsFilters(x, y, chi_cutoff);
               }
               catch (std::exception& e)
               {
                  std::cerr << "kmer " + k.sequence() + " failed chisq test with error: " + e.what();
                  pass = 1;
               }

               if (pass)
               {
#ifdef PANGWAS_DEBUG
                  std::cerr << "kmer " + k.sequence() + " seems significant\n";
#endif
                  k.add_x(x);
                  kmer_lines.push_back(k);
               }
            }
         }
      }

#ifdef NO_THREAD
      if (kmer_lines.size() == 1)
      {
         logisticTest(kmer_lines[0], y);

         if (kmer_lines[0].p_val() < log_cutoff)
         {
            std::cout << kmer_lines[0];
         }
      }
#else
      // Thread from here...
      // TODO if slow, can thread with multiple tests per thread
      std::vector<std::thread> threads;
      threads.reserve(kmer_lines.size());

      for (unsigned int i = 0; i<kmer_lines.size(); ++i)
      {
         // Association test
         // Note threads must be passed values as they are copied
         // std::reference_wrapper allows references to be passed
         threads.push_back(std::thread(logisticTest, std::ref(kmer_lines[i]), std::cref(y)));
      }

      for (unsigned int i = 0; i<threads.size(); ++i)
      {
         // Rejoin in order
         threads[i].join();

         // Print in order when all threads complete
         if (kmer_lines[i].p_val() < log_cutoff)
         {
            std::cout << kmer_lines[i];
         }
      }
      // ...to here
#endif
   }

   std::cerr << "Done.\n";

}

