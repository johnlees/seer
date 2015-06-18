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
   int continuous_phenotype = continuousPhenotype(samples);

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

   void (*mdsLinearFunc)(Kmer&, const arma::vec&, const arma::mat&) = &linearTest;
   void (*linearFunc)(Kmer&, const arma::vec&) = &linearTest;
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
            else if (passFilters(parameters, k, samples, y, continuous_phenotype))
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
            if (continuous_phenotype)
            {
               linearTest(kmer_lines[0], y, mds);
            }
            else
            {
               logisticTest(kmer_lines[0], y, nr_opt, mds);
            }
         }
         else
         {
            if (continuous_phenotype)
            {
               linearTest(kmer_lines[0], y);
            }
            else
            {
               logisticTest(kmer_lines[0], y, nr_opt);
            }
         }

         if (kmer_lines[0].p_val() < parameters.log_cutoff)
         {
            std::cout << kmer_lines[0];
            if (parameters.print_samples)
            {
               std::vector<std::string> samples_found = kmer_lines[0].occurrence_vector();

               // Doing this for all samples leaves trailing whitespace, so
               // write the last sample separately
               if (samples_found.size() > 1)
               {
                  std::copy(samples_found.begin(), samples_found.end() - 1, std::ostream_iterator<std::string>(std::cout, "\t"));
               }
               std::cout << *samples_found.end();
               std::cout << "\n";
            }
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
            if (continuous_phenotype)
            {
               threads.push_back(std::thread(mdsLinearFunc, std::ref(kmer_lines[i]), std::cref(y), std::cref(mds)));
            }
            else
            {
               threads.push_back(std::thread(mdsLogitFunc, std::ref(kmer_lines[i]), std::cref(y), nr_opt, std::cref(mds)));
            }
         }
         else
         {
            if (continuous_phenotype)
            {
               threads.push_back(std::thread(linearFunc, std::ref(kmer_lines[i]), std::cref(y)));
            }
            else
            {
               threads.push_back(std::thread(logitFunc, std::ref(kmer_lines[i]), std::cref(y), nr_opt));
            }
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
            if (parameters.print_samples)
            {
               std::vector<std::string> samples_found = kmer_lines[i].occurrence_vector();

               // Doing this for all samples leaves trailing whitespace, so
               // write the last sample separately
               if (samples_found.size() > 1)
               {
                  std::copy(samples_found.begin(), samples_found.end() - 1, std::ostream_iterator<std::string>(std::cout, "\t"));
               }
               std::cout << samples_found.back();
               std::cout << "\n";
            }

         }
      }
      // ...to here
#endif
   }

   std::cerr << "Done.\n";

}

