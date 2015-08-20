/*
 * Header file for combineKmers.cpp
 * Takes union of dsk counted kmers
 *
 */

// C++ stl includes
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iterator>

// Library includes
#include <H5Cpp.h>

// Global defaults
const int default_kmer_size = 31;
const std::string dsk_group_name = "dsk/solid";

// Function prototypes
std::string readH5StrAttr(const H5::H5File& h5_file, const std::string group_name, const std::string attr_name);
