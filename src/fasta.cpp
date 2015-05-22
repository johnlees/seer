/*
 * File: fasta.cpp
 *
 * Helper functions for the fasta sequence class
 *
 */

#include "fasta.hpp"

// Read in sequences to memory from file
Fasta::Fasta(const std::string& obj_name, const std::string& filename)
   :name(obj_name)
{
   std::ifstream ist(filename.c_str());

   if (!ist)
   {
      throw std::runtime_error("Could not open fasta file " + filename + "\n");
   }

   std::string line_in;
   std::string seq, contig_name = "";

   std::getline(ist, line_in);

   if (ist.eof())
   {
       throw std::runtime_error("Could not read fasta file " + filename + "\n");
   }
   // Read header line, which starts with >
   else if (line_in[0] == '>')
   {
      contig_name = line_in.substr(1);
   }
   else
   {
      throw std::runtime_error("Error reading fasta file " + filename + "\n"
         "Header should start with '>', but has:\n" + line_in + "\n");
   }

   while (!ist.eof())
   {
      // New contig, new sequence
      if (ist.peek() == '>')
      {
         sequence_names.push_back(contig_name);
         sequences.push_back(seq);

         seq = "";
         std::getline(ist, line_in);
         contig_name = line_in.substr(1);
      }
      else
      {
         // Concatenate sequence lines
         std::getline(ist, line_in);
         seq += line_in;
      }
   }

   // Add in final contig
   sequence_names.push_back(contig_name);
   sequences.push_back(seq);

}

// All exact matches to search sequence in all sequences in fasta
std::vector<Mapping> Fasta::hasSeq(const std::string& search)
{
   std::vector<Mapping> results;
   Mapping hit;

   // Search each sequence in fasta
   std::vector<std::string>::iterator name_it = Fasta::sequence_names.begin();
   for (std::vector<std::string>::iterator it = Fasta::sequences.begin(); it != Fasta::sequences.end(); ++it)
   {
      size_t hit_position = it->find(search);

      // Search allowing for multiple hits
      while (hit_position != std::string::npos)
      {
         hit.sequence_name = *name_it;
         hit.position = hit_position;
         results.push_back(hit);

         hit_position = it->find(search, hit_position + 1);
      }

      // Move forward in name vector
      ++name_it;
   }

   return results;
}

// Prints results to hasSeq
// a nicer output than hasSeq: use this function
void Fasta::printMappings(std::ostream &os, const std::string& search)
{
   std::vector<Mapping> hits = hasSeq(search);

   if (hits.size() > 0)
   {
      for (std::vector<Mapping>::iterator it = hits.begin(); it != hits.end(); ++it)
      {
         os << "\t" << Fasta::name << ":" << it->sequence_name << ":" << it->position;
      }
   }
}

