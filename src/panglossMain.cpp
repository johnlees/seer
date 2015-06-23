/*
 * File: panglossMain.cpp
 *
 * Control process of pangloss (mds of kmers)
 *
 */

#include "pangloss.hpp"

// globals
std::default_random_engine rand_gen;

// $1 file location, $2 file name, $3 file ending
const std::regex file_format_e ("^(.+)\\/(.+)\\.([^\\.]+)$");
const std::regex file_format_within_e ("^(.+)\\.([^\\.]+)$"); // When in same directory, but not with a ./ prefix

int main (int argc, char *argv[])
{
   // STRUCTURE
   // Read in command line
   //    dsm and pheno file
   //    Matrix size
   //    No. of prin components
   //    Have a no_filter option
   //    Output prefix, with default (same as input)
   //
   // Read a dsm line, convert to stringstream
   //    Parse as kmer
   //    Write out if passes filters
   //
   // Build large matrix using reservoir sampler
   //    Reserve enough space
   //    While under reserved space, add kmer to vector
   //    Then reservoir sample
   //
   // Convert kmers into matrix
   //
   // PCoA on matrix
   //    Save matrix of top 3 eigenvals (parameter input)

   // Program description
   std::cerr << "pangloss: pan-genome limitation of structure sensitivity\n";

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   boost::program_options::variables_map vm;
   if (argc == 1)
   {
      std::cerr << "Usage: pangloss -k dsm.txt.gz -p data.pheno\n\n"
         << "For full option details run pangloss -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, vm))
   {
      return 1;
   }

   // Normal operation
   if (!vm.count("mds_concat"))
   {
      // Check compulsory paramters provided
      if (!vm.count("kmers") || !vm.count("pheno"))
      {
         std::cerr << "Either -k and -p must be provided, or --mds_concat\n";
         return 1;
      }

      // Open .pheno file, parse into vector of samples
      std::vector<Sample> samples;
      readPheno(vm["pheno"].as<std::string>(), samples);
      arma::vec y = constructVecY(samples);
      int continuous_phenotype = continuousPhenotype(samples);

      cmdOptions parameters = verifyCommandLine(vm, samples);

      // Open the dsm kmer ifstream, and read through the whole thing
      igzstream kmer_file;
      openDsmFile(kmer_file, parameters.kmers.c_str());

      // Set up output files
      ogzstream filtered_file;
      std::string output_file_name, dsm_file_name;
      if (parameters.filter)
      {
         if (vm.count("output"))
         {
            output_file_name = parameters.output + ".kmers.gz";
         }
         else
         {
            output_file_name = std::regex_replace(parameters.kmers, file_format_e, std::string("$1/filtered.$2.gz"));
            if (output_file_name == parameters.kmers) // If first match fails
            {
               output_file_name = std::regex_replace(parameters.kmers, file_format_within_e, std::string("filtered.$1.gz"));
            }

         }

         filtered_file.open(output_file_name.c_str());
      }

      if (vm.count("output"))
      {
         dsm_file_name = parameters.output + ".dsm";
      }
      else
      {
         dsm_file_name = std::regex_replace(parameters.kmers, file_format_e, std::string("$1/$2.dsm"));
         if (dsm_file_name == parameters.kmers) // If first match fails
         {
            dsm_file_name = std::regex_replace(parameters.kmers, file_format_within_e, std::string("$1.dsm"));
         }
      }

      // vector of subsampled kmers
      std::vector<arma::vec> dsm_kmers;
      dsm_kmers.reserve(parameters.size);

      long int kmer_index = 0;
      while (kmer_file)
      {
         // Allow output of entire dsm line
         std::string dsm_line;
         std::getline(kmer_file, dsm_line);
         std::stringstream dsm_stream(dsm_line);

         Kmer k;
         if (kmer_file)
         {
            dsm_stream >> k;

            // apply filters here
            if (!parameters.filter)
            {
               k.add_x(constructVecX(k, samples));
            }
            else if (passFilters(parameters, k, samples, y, continuous_phenotype))
            {
#ifdef PANGWAS_DEBUG
               std::cerr << "kmer " + k.sequence() + " seems significant\n";
#endif
               filtered_file << dsm_line << "\n";
            }

            // kmer has passed basic filters, so is a candidate for mds
            // subsampling
            if (k.has_x())
            {
               // Resevoir sampler for parameters.size kmers
               kmer_index++;
               if (dsm_kmers.size() < dsm_kmers.capacity())
               {
                  dsm_kmers.push_back(k.get_x());
               }
               else
               {
                  std::uniform_int_distribution<int> dist(0, kmer_index);
                  int r = dist(rand_gen);
                  if (r < parameters.size)
                  {
                     dsm_kmers[r] = k.get_x();
                  }
               }
            }
         }
      }

      // Convert into arma::mat before MDS
      arma::mat subsampledMatrix(samples.size(), dsm_kmers.size());
      for (unsigned int i = 0; i < dsm_kmers.size(); ++i)
      {
         subsampledMatrix.col(i) = dsm_kmers[i];
      }

      // Write output
      if (vm.count("no_mds"))
      {
         writeMDS(dsm_file_name, subsampledMatrix);
      }
      else
      {
         // Run metric MDS, then output to file
         writeMDS(dsm_file_name, metricMDS(subsampledMatrix, parameters.pc, parameters.num_threads));
      }

      std::cerr << "Done.\n";
      if (!vm.count("no_mds"))
      {
         if (parameters.filter)
         {
            std::cerr << "You may now want to run:\n\tpangwas -k " << output_file_name
               << " -p " << vm["pheno"].as<std::string>() << " --struct " << dsm_file_name
               << " > significant_kmers.txt\n";
         }
         else
         {
            std::cerr << "You may now want to run:\n\tpangwas -k " << parameters.kmers
               << " -p " << vm["pheno"].as<std::string>() << " --struct " << dsm_file_name
               << " > significant_kmers.txt\n";
         }
      }
   }
   // Matrices already subsampled. Concatenate and perform MDS
   else
   {
      std::string matrix_input = vm["mds_concat"].as<std::string>();

      // Check number of principal coordinates ok
      cmdOptions mdsOptions;
      verifyMDSOptions(mdsOptions, vm);

      std::string dsm_file_name;
      if (fileStat(matrix_input))
      {
         if (vm.count("output"))
         {
            dsm_file_name = vm["output"].as<std::string>();
         }
         else
         {
            dsm_file_name = matrix_input + ".dsm";
         }
         // Run metric MDS, then output to file
         writeMDS(dsm_file_name, metricMDS(readMDSList(matrix_input), mdsOptions.pc, mdsOptions.num_threads));
      }
      else
      {
         throw std::runtime_error("Couldn't open list of matrices " + matrix_input);
      }

      std::cerr << "Done.\n";

      std::cerr << "Output written to " << dsm_file_name << "\n"
         << "Use this as the --struct option of pangwas\n";
   }
}

