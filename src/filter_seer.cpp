/*
 * filter_seer.cpp
 * Filters significant kmers produced by seer
 *
 */

#include "filter_seer.hpp"

int main (int argc, char *argv[])
{
   // Program description
   std::cerr << "filter_seer: post filtering of significant kmers\n";

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   boost::program_options::variables_map vm;
   if (argc == 1)
   {
      std::cerr << "Usage: filter_seer -k seer_kmers.txt --pos_beta > filtered_output.txt\n\n"
         << "For full option details run 'filter_seer -h'\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, vm))
   {
      return 1;
   }

   cmdOptions options = processCmdLine(vm);

   // Define sorting function object
   sortSigKmer sort_kmers;
   if (options.sort_field != "")
   {
      sort_kmers.addSortField(options.sort_field);
   }

   sortSigKmer substr_sort;
   if (!options.substr_kmers)
   {
      substr_sort.addSortField("sequence");
   }

   std::list<Significant_kmer> filtered_kmers;

   // Open input file, check opened ok
   try
   {
      std::ifstream kmers_in(options.input_file.c_str());
      if (!kmers_in)
      {
         throw std::runtime_error("Could not open input file " + options.input_file);
      }
      else
      {
         while (kmers_in)
         {
            // Read in each kmer
            Significant_kmer kmer;
            kmers_in >> kmer;

            if (kmers_in)
            {
               // Filter fields. MAF needs to include max, beta needs modulus
               if (options.maf_filter && (kmer.maf() < options.maf_filter || 1-kmer.maf() < options.maf_filter))
               {
                  continue;
               }
               else if (!options.neg_beta && kmer.beta() < 0)
               {
                  continue;
               }
               else if (options.beta_filter && fabs(kmer.beta()) < options.beta_filter)
               {
                  continue;
               }
               else if (options.chi_filter && kmer.unadj() > options.chi_filter)
               {
                  continue;
               }
               else if (options.p_filter && kmer.p_val() > options.p_filter)
               {
                  continue;
               }
               // If sorted, store in mem. Otherwise print immediately
               else if (!sort_kmers && !substr_sort)
               {
                  std::cout << kmer << std::endl;
               }
               else
               {
                  filtered_kmers.push_back(kmer);
               }
            }
         }
         // Remove substr if necessary
         if (substr_sort)
         {
            // Sort largest to smallest
            filtered_kmers.sort(substr_sort);

            auto list_it = filtered_kmers.begin();
            while (list_it != filtered_kmers.end())
            {
               // Check in all the larger kmers that this one isn't
               // a substring
               bool found = 0;
               for (auto sub_it = filtered_kmers.begin(); sub_it != list_it; ++sub_it)
               {
                  if (sub_it->sequence().find(list_it->sequence()) != std::string::npos)
                  {
                     list_it = filtered_kmers.erase(list_it);
                     found = 1;
                     break;
                  }
               }

               if (!found)
               {
                  ++list_it;
               }
            }
         }

         // If sorting, sort, then print
         if (sort_kmers)
         {
            filtered_kmers.sort(sort_kmers);
         }

         if (sort_kmers || substr_sort)
         {
            for (auto out_it = filtered_kmers.begin(); out_it != filtered_kmers.end(); ++out_it)
            {
               std::cout << *out_it << std::endl;
            }
         }

         std::cerr << "Done." << std::endl;
      }
   }
   catch (std::exception& e)
   {
      std::cerr << "Error: " << e.what() << std::endl;
   }

   return(0);
}

