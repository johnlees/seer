/*
 * File: seerMain.cpp
 *
 * Reads command line input to seer, and controls relevant functions
 *
 */

#include "seer.hpp"

int main (int argc, char *argv[])
{
   // Program description
   std::cerr << "seer: sequence element enrichment analysis\n";

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   boost::program_options::variables_map vm;
   if (argc == 1)
   {
      std::cerr << "Usage: seer -k dsm.txt.gz -p data.pheno --struct mds.dsm.txt.gz\n\n"
         << "For full option details run seer -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, vm))
   {
      return 1;
   }

   // Open .pheno file, parse into vector of samples
   std::vector<Sample> samples;
   std::unordered_map<std::string,int> sample_map;
   readPheno(vm["pheno"].as<std::string>(), samples, sample_map);
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

   // Set up covariates
   if (vm.count("covar_file") && vm.count("covar_list"))
   {
      // File reading/parsing may fail
      try
      {
         arma::mat covariate_matrix =
            parseCovars(vm["covar_file"].as<std::string>(), vm["covar_list"].as<std::string>());

         if (covariate_matrix.n_rows != samples.size())
         {
            throw std::runtime_error("Covariate row size does not match number of samples");
         }

         if (use_mds)
         {
            mds = arma::join_rows(mds, covariate_matrix);
         }
         else
         {
            mds = covariate_matrix;
            use_mds = 1;
         }
      }
      catch (std::exception& e)
      {
         std::cerr << "Could not process covariates: " << std::endl;
         std::cerr << e.what() << std::endl;
      }
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

   // Disambiguate overloaded logistic functions by the type of parameter they
   // take
   void (*mdsLogitFunc)(Kmer&, const arma::vec&, const unsigned int nr, const arma::mat&) = &logisticTest;
   void (*logitFunc)(Kmer&, const arma::vec&, const unsigned int nr) = &logisticTest;

   void (*mdsLinearFunc)(Kmer&, const arma::vec&, const arma::mat&) = &linearTest;
   void (*linearFunc)(Kmer&, const arma::vec&) = &linearTest;

   // Error check command line options
   cmdOptions parameters = verifyCommandLine(vm, samples);

   // Open the dsm kmer ifstream, and read through the whole thing
   igzstream kmer_file;
   openDsmFile(kmer_file, parameters.kmers);

   while (kmer_file)
   {
      // Parse a set of dsm lines
      std::vector<Kmer> kmer_lines;
      kmer_lines.reserve(parameters.num_threads);

      Kmer k;
      while(kmer_lines.size() < parameters.num_threads && kmer_file)
      {
         kmer_file >> k;
         if (kmer_file)
         {
            k.add_x(sample_map, samples.size());

            // apply filters here
            if (!parameters.filter)
            {
               kmer_lines.push_back(k);
            }
            else if (passFilters(parameters, k, samples, y, continuous_phenotype))
            {
#ifdef SEER_DEBUG
               std::cerr << "kmer " + k.sequence() + " seems significant\n";
#endif
               kmer_lines.push_back(k);
            }
         }
      }

      // Thread from here...
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
            // Caclculate chisq value if not already done so in filtering
            if (kmer_lines[i].chi_p_val() == kmer_chi_pvalue_default)
            {
               if (continuous_phenotype)
               {
                  kmer_lines[i].chi_p_val(welchTwoSamplet(kmer_lines[i].get_x(), y));
               }
               else
               {
                  kmer_lines[i].chi_p_val(chiTest(kmer_lines[i].get_x(), y));
               }
            }

            std::cout << kmer_lines[i];
            if (parameters.print_samples)
            {
               std::vector<std::string> samples_found = kmer_lines[i].occurrence_vector();
               std::cout << "\t";
               // Doing this for all samples leaves trailing whitespace, so
               // write the last sample separately
               if (samples_found.size() > 1)
               {
                  std::copy(samples_found.begin(), samples_found.end() - 1, std::ostream_iterator<std::string>(std::cout, "\t"));
               }
               std::cout << samples_found.back();
            }
            std::cout << "\n";

         }
      }
      // ...to here
   }

   std::cerr << "Done.\n";
}

