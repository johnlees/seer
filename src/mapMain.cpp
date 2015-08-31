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
      std::cerr << "Usage: map_back -k seer_output.txt -r references.txt > mappings.txt\n\n"
         << "For full option details run map_back -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, vm))
   {
      return 1;
   }

   // Set number of threads
   size_t num_threads = vm["threads"].as<size_t>();
   if (num_threads < 1)
   {
      num_threads = 1;
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
      // Set up list of threads, and mutex lock on std out
      std::mutex out_mtx;
      std::list<std::future<void>> threads;

      Significant_kmer sig_kmer;
      while(kmer_file)
      {
         kmer_file >> sig_kmer;

         // Check the read into sig_kmer hasn't reached end of file
         if (!kmer_file.eof())
         {
            assert(num_threads >= 1);

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
                  // Thread each search (i.e. per sample per kmer)
                  if (num_threads > 1)
                  {
                     // Check if four searches are running. If so, wait for the
                     // next one to finish. Wait (for a max of 100ms) so this
                     // thread doesn't consume processing
                     waitForThreads(threads, num_threads - 1);

                     // Start a new thread asynchronously, at the back of the
                     // queue
                     threads.push_back(std::async(std::launch::async,
                              &Fasta::printMappings, *all_names_it, std::ref(std::cout), sig_kmer.sequence(), std::ref(out_mtx)));
                  }
                  else
                  {
                     all_names_it->printMappings(std::cout, sig_kmer.sequence(), out_mtx);
                  }

                  ++search_names_it;
                  if (search_names_it == search_names.end())
                  {
                     break;
                  }
               }
            }

            // Tab between every sample, line break after every kmer
            if (num_threads > 1)
            {
               waitForThreads(threads, 0);
            }

            std::cout << std::endl;
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

// Waits until there are #spaces threads left unfinished
void waitForThreads(std::list<std::future<void>>& thread_list, const size_t leave_running)
{
   auto list_it = thread_list.begin();
   while (thread_list.size() > leave_running)
   {
      auto status = list_it->wait_for(std::chrono::milliseconds(thread_wait));
      if (status == std::future_status::ready)
      {
         list_it = thread_list.erase(list_it);
      }
      else if (list_it++ == thread_list.end())
      {
         list_it = thread_list.begin();
      }
   }
}
