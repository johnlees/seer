/*
 * File: pangwasMain.cpp
 *
 * Reads command line input to pangwas, and controls relevant functions
 *
 */

#include "pangwas.h"

namespace po = boost::program_options; // Save some typing

int main (int argc, char *argv[])
{
   po::variables_map vm;

   // Do parsing and checking of command line params
   // If no input options, give quick usage rather than full help
   if (argc == 1)
   {
      std::cerr << "Usage: pangwas -k dsm.txt -p data.pheno -o output\n\n"
         << "For full option details run pangwas -h\n";
      return 0;
   }
   else if (parseCommandLine(argc, argv, &vm))
   {
      return 1;
   }

   if (vm.count("compression"))
   {
      std::cout << "Compression level was set to " << vm["compression"].as<int>() << ".\n";
   }
   else
   {
      std::cout << "Compression level was not set.\n";
   }
}

// Use boost::program_options to parse command line input
int parseCommandLine (int argc, char *argv[], po::variables_map *vm)
{
   int failed = 0;

   po::options_description desc("Options");
   desc.add_options()
    ("help,h", "produce help message")
    ("compression", po::value<int>(), "set compression level");

   try
   {
      po::store(po::parse_command_line(argc, argv, desc), *vm);
      po::notify(*vm);
   }
   catch (po::error& e)
   {
      std::cerr << "Error in command line input: " << e.what() << "\n\n";
      std::cerr << desc << "\n";
      return 1;
   }

   // Print long help message
   if (vm->count("help"))
   {
      std::cerr << desc << "\n";
      failed = 1;
   }
   else
   {
      failed = 0;
   }

   return failed;
}


