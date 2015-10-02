/*
 * File: seerIO.cpp
 *
 * Controls: reading in pheno to sample object
 * Open (zipped) dsm file
 * Internal format conversion
 *
 * NB: some input and output used overridden operators in
 * <class>.cpp files
 *
 */
#include "seercommon.hpp"

std::regex gzipped(".+\\.gz"); // matches .gz at the end of a string (as regex match
                               // matches a whole string)

/* .pheno files
 * space separated:
 * FID, IID, phenotype (1 control, 2 case)
 */
void readPheno(const std::string& filename, std::vector<Sample>& samples, std::unordered_map<std::string,int>& sample_map)
{
   std::ifstream ist(filename.c_str());

   if (!ist)
   {
      throw std::runtime_error("Could not open pheno file " + filename + "\n");
   }

   Sample s;
   int sample_index = 0;
   while (ist >> s)
   {
      samples.push_back(s);

      sample_map[s.iid()] = sample_index;
      sample_index++;
   }
}

// Open dsm files, which are possibly zipped
void openDsmFile(igzstream& dsm_stream, const std::string& file_name)
{
   // Check for a .gz extension
   if (!std::regex_match(file_name, gzipped))
   {
      // Warn
      std::cerr << "WARNING: Input file " + file_name
         + " is not gzip compressed, which is recommended\n";
   }

   dsm_stream.open(file_name.c_str()); // Push binary file into buffer

   // Set stream buffer of istream to the one just opened, and check ok
   if (!dsm_stream.good())
   {
      throw std::runtime_error("Could not open kmer file " + file_name + "\n");
   }

}

arma::vec constructVecY(const std::vector<Sample>& samples)
{
   arma::vec y;
   y.zeros(samples.size());

   for (unsigned int i = 0; i < samples.size(); ++i)
   {
      y[i] = samples[i].pheno();
   }

   return y;
}

void writeMDS(const std::string& file_name, const arma::mat& MDS)
{
   MDS.save(file_name, arma::hdf5_binary);
}

void writeDistances(const std::string& file_name, const arma::mat& distances)
{
   distances.save(file_name, arma::csv_ascii);
}

arma::mat readMDS(const std::string& file_name)
{
   arma::mat MDS;

   if (fileStat(file_name))
   {
      MDS.load(file_name);
   }
   else
   {
      throw std::runtime_error("Problem with loading MDS input file " + file_name);
   }

   return MDS;
}

arma::mat readMDSList(const std::string& filename)
{
   std::ifstream ist(filename.c_str());
   if (!ist)
   {
      throw std::runtime_error("Could not open mds list file " + filename + "\n");
   }
   else
   {
      std::cerr << "Reading subsampled matrices from " + filename + "\n";
   }

   std::string matrix_file;
   arma::mat combined_matrix;
   int i = 0;
   while (ist >> matrix_file)
   {
      combined_matrix = join_rows(combined_matrix, readMDS(matrix_file));

      std::cerr << "Joined matrix " << ++i << "\n";
   }

   return combined_matrix;
}

// Check for file existence
int fileStat(const std::string& filename)
{
   struct stat buffer;
   int success = 1;

   if (stat (filename.c_str(), &buffer) != 0)
   {
      std::cerr << "Can't stat input file: " << filename << "\n";

      success = 0;
   }

   return success;
}

