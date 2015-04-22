/*
 * File: pangwasMain.cpp
 *
 * Reads command line input to pangwas, and controls relevant functions
 *
 */

#include "pangwas.hpp"

int main (int argc, char *argv[])
{
   // Program description
   std::cerr << "pangwas: pan-genome wide association study using kmers\n";

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   boost::program_options::variables_map vm;
   if (argc == 1)
   {
      std::cerr << "Usage: pangwas -k dsm.txt.gz -p data.pheno --struct mds.dsm.txt.gz\n\n"
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
   arma::vec y = constructVecY(samples);

   // Get mds values
   arma::mat mds;
   int use_mds = 0;
   if (vm.count("struct"))
   {
      mds = readMDS(vm["struct"].as<std::string>());
      use_mds = 1;

   }

   unsigned int nr_opt;
   if (vm.count("optimiser"))
   {
      if (vm["optimiser"].as<std::string>() == "newton-raphson")
      {
         nr_opt = 1;
      }
      else
      {
         nr_opt = 0;
      }
   }

#ifndef NO_THREAD
   // Disambiguate overloaded logistic functions by the type of parameter they
   // take
   void (*mdsLogitFunc)(Kmer&, const arma::vec&, const unsigned int nr, const arma::mat&) = &logisticTest;
   void (*logitFunc)(Kmer&, const arma::vec&, const unsigned int nr) = &logisticTest;
#endif

   // Error check command line options
   cmdOptions parameters = verifyCommandLine(vm, samples);

   // Open the dsm kmer ifstream, and read through the whole thing
   igzstream kmer_file;
   openDsmFile(kmer_file, parameters.kmers);

   while (kmer_file)
   {
      // Parse a set of dsm lines
      // TODO this could also be threaded?
      std::vector<Kmer> kmer_lines;
      kmer_lines.reserve(parameters.num_threads);

      Kmer k;
      while(kmer_lines.size() < parameters.num_threads && kmer_file)
      {
         if (kmer_file)
         {
            kmer_file >> k;

            // apply filters here
            if (!parameters.filter)
            {
               k.add_x(constructVecX(k, samples));
               kmer_lines.push_back(k);
            }
            else if (passFilters(parameters, k, samples, y))
            {
#ifdef PANGWAS_DEBUG
               std::cerr << "kmer " + k.sequence() + " seems significant\n";
#endif
               kmer_lines.push_back(k);
            }
         }
      }

#ifdef NO_THREAD
      if (kmer_lines.size() == 1)
      {
         if (use_mds)
         {
            logisticTest(kmer_lines[0], y, nr_opt, mds);
         }
         else
         {
            logisticTest(kmer_lines[0], y, nr_opt);
         }

         if (kmer_lines[0].p_val() < parameters.log_cutoff)
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
         if (use_mds)
         {
            threads.push_back(std::thread(mdsLogitFunc, std::ref(kmer_lines[i]), std::cref(y), nr_opt, std::cref(mds)));
         }
         else
         {
            threads.push_back(std::thread(logitFunc, std::ref(kmer_lines[i]), std::cref(y), nr_opt));
         }
      }

      for (unsigned int i = 0; i<threads.size(); ++i)
      {
         // Rejoin in order
         threads[i].join();

         // Print in order when all threads complete
         if (kmer_lines[i].p_val() < parameters.log_cutoff)
         {
            std::cout << kmer_lines[i];
         }
      }
      // ...to here
#endif
   }

   std::cerr << "Done.\n";

}

