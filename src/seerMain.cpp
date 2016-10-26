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

   if (vm.count("pheno"))
   {
      readPheno(vm["pheno"].as<std::string>(), samples, sample_map);
   }
   else
   {
      throw std::runtime_error("--pheno option is compulsory");
   }

   arma::vec y = constructVecY(samples);
   int continuous_phenotype = continuousPhenotype(samples);

   // Get mds values
   arma::mat mds;
   int use_mds = 0;
   if (vm.count("struct"))
   {
      mds = readMDS(vm["struct"].as<std::string>(), samples);
      use_mds = 1;

      if (mds.n_rows != samples.size())
      {
         throw std::runtime_error("Number of rows in MDS matrix does not match number of samples");
      }
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

   // Disambiguate overloaded logistic functions by the type of parameter they
   // take
   void (*mdsLogitFunc)(Kmer&, const arma::vec&, const double, const arma::mat&) = &logisticTest;
   void (*logitFunc)(Kmer&, const arma::vec&, const double) = &logisticTest;

   void (*mdsLinearFunc)(Kmer&, const arma::vec&, const double, const arma::mat&) = &linearTest;
   void (*linearFunc)(Kmer&, const arma::vec&, const double) = &linearTest;

   // Error check command line options
   cmdOptions parameters = verifyCommandLine(vm, samples);

   // calculate the null log-likelihood
   Kmer null_kmer;
   arma::mat x(y.n_rows, 1, arma::fill::ones);
   if (use_mds)
   {
      x = join_rows(x, mds);
   }

   double null_ll = nullLogLikelihood(x, y, continuous_phenotype);

   // Open the dsm kmer ifstream, and read through the whole thing
   igzstream kmer_file;
   openDsmFile(kmer_file, parameters.kmers);

   // Write a header
   std::string header = "sequence\tmaf\tchisq_p_val\twald_p_val\tlrt_p_val\tbeta\tse";
   if (use_mds)
   {
      std::cout << header;
      for (unsigned int i = 1; i <= mds.n_cols; ++i)
      {
         std::cout << "\tcovar" << i << "_p";
      }
      std::cout << "\tcomments";
   }
   else
   {
      std::cout << header << "\tcomments";
   }
   if (parameters.print_samples)
   {
      std::cout << "\tsamples_present";
   }
   std::cout << std::endl;

   long int input_line = 0;
   long int tested_kmers = 0;
   long int significant_kmers = 0;
   while (kmer_file)
   {
      // Parse a set of dsm lines
      std::vector<Kmer> kmer_lines;
      kmer_lines.reserve(parameters.num_threads);

      Kmer k;
      while(kmer_lines.size() < parameters.num_threads && kmer_file)
      {
         kmer_file >> k;
         k.set_line_nr(++input_line);

         if (kmer_file)
         {
            k.add_x(sample_map, samples.size());

            // apply filters here
            if (!parameters.filter)
            {
               kmer_lines.push_back(k);
               tested_kmers++;
            }
            else if (passBasicFilters(parameters, k) && passStatsFilters(parameters, k, y, continuous_phenotype))
            {
#ifdef SEER_DEBUG
               std::cerr << "kmer " + k.sequence() + " seems significant\n";
#endif
               kmer_lines.push_back(k);
               tested_kmers++;
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
               threads.push_back(std::thread(mdsLinearFunc, std::ref(kmer_lines[i]), std::cref(y), null_ll, std::cref(mds)));
            }
            else
            {
               threads.push_back(std::thread(mdsLogitFunc, std::ref(kmer_lines[i]), std::cref(y), null_ll, std::cref(mds)));
            }
         }
         else
         {
            if (continuous_phenotype)
            {
               threads.push_back(std::thread(linearFunc, std::ref(kmer_lines[i]), std::cref(y), null_ll));
            }
            else
            {
               threads.push_back(std::thread(logitFunc, std::ref(kmer_lines[i]), std::cref(y), null_ll));
            }
         }
      }

      for (unsigned int i = 0; i<threads.size(); ++i)
      {
         // Rejoin in order
         threads[i].join();

         // Print in order when all threads complete
         if (kmer_lines[i].p_val() < parameters.log_cutoff || kmer_lines[i].lrt_p_val() < parameters.log_cutoff)
         {
            significant_kmers++;
            // Caclculate chisq value if not already done so in filtering
            if (kmer_lines[i].unadj() == kmer_chi_pvalue_default)
            {
               if (continuous_phenotype)
               {
                  kmer_lines[i].unadj_p_val(welchTwoSamplet(kmer_lines[i], y));
               }
               else
               {
                  kmer_lines[i].unadj_p_val(chiTest(kmer_lines[i], y));
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
            std::cout << std::endl;

         }
      }
      // ...to here
   }

   std::cerr << "Read " << input_line - 1 << " total k-mers. Of these:\n";
   std::cerr << "\tPre-filtered " << input_line - tested_kmers - 1 << " k-mers\n";
   std::cerr << "\tTested " << tested_kmers << " k-mers\n";
   std::cerr << "\tPrinted " << significant_kmers << " k-mers\n";
   std::cerr << "Done.\n";
}

