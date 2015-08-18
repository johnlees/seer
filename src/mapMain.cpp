/*
 * File: mapMain.cpp
 *
 * Reports exact matches of kmers to assemblies
 *
 */

#include "map_back.hpp"

int main (int argc, char *argv[])
{
   // Program description
   std::cerr << "map_back: context of significant kmers\n";

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   boost::program_options::variables_map vm;
   if (argc == 1)
   {
      std::cerr << "Usage: map_back -k pangwas_output.txt -r references.txt > mappings.txt\n\n"
         << "For full option details run map_back -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, vm))
   {
      return 1;
   }

   // Read all sequences into memory as Fasta objects
   std::cerr << "Reading reference sequences into memory...\n";
   std::vector<Fasta> sequence_cache = readSequences(vm["references"].as<std::string>());

   // Loop through significant kmers
   std::cerr << "Now mapping significant kmers...\n";
   std::ifstream kmer_file(vm["kmers"].as<std::string>().c_str());
   if (!kmer_file)
   {
      throw std::runtime_error("Could not open kmer_file " + vm["kmers"].as<std::string>() + "\n");
   }
   else
   {
      Significant_kmer sig_kmer;
      while(kmer_file)
      {
         kmer_file >> sig_kmer;

         // Check the read into sig_kmer hasn't reached end of file
         if (!kmer_file.eof())
         {
            std::cout << sig_kmer.sequence();

            // sig_kmer samples and sample cache are sorted in the same order, so
            // can go through linearly
            std::vector<std::string> search_names = sig_kmer.samples_found();
            std::vector<std::string>::iterator search_names_it = search_names.begin();

            for (std::vector<Fasta>::iterator all_names_it = sequence_cache.begin(); all_names_it != sequence_cache.end(); ++all_names_it)
            {
               // For each sample we know the kmer is in, print all matches to
               // the kmer
               if (all_names_it->get_name() == *search_names_it)
               {
                  all_names_it->printMappings(std::cout, sig_kmer.sequence());

                  ++search_names_it;
                  if (search_names_it == search_names.end())
                  {
                     break;
                  }
               }
            }
            // Tab between every sample, line break after every kmer
            std::cout << "\n";
         }
      }

      std::cerr << "Done.\n";
   }

}

// Stores all fasta sequences in a vector of Fasta objects
std::vector<Fasta> readSequences(const std::string& reference_file)
{
   std::vector<Fasta> sequences;

   std::ifstream ifs(reference_file.c_str());
   if (!ifs)
   {
      throw std::runtime_error("Could not open reference file: " + reference_file + "\n");
   }
   else
   {
      std::string name, fasta_file = "";
      while(ifs)
      {
         ifs >> name >> fasta_file;
         sequences.push_back(Fasta(name, fasta_file));
      }
   }

   std::sort(sequences.begin(), sequences.end(), Fasta::compareFasta);
   return sequences;
}
