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
std::regex covar_regex("^(\\d+)$");
std::regex q_covar_regex("^(\\d+)q$");

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

   // Always keep samples sorted, for consistency between programs
   std::sort(samples.begin(), samples.end(), Sample::compareSamples);
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

// Opens a text file containing covariates, and reads the specified columns
// into a matrix. Categorical covariates use dummy coding
arma::mat parseCovars(const std::string& file, const std::string& columns)
{
   // Check file exists
   if (!fileStat(file))
   {
      throw std::runtime_error("Covariate file " + file + " not found");
   }

   // First parse columns required
   std::vector<std::tuple<int,bool>> column_list = parseCovarColumns(columns);

   // Open file, and read in required columns
   std::vector<Covar> covars_in(column_list.size()); // Vector of covariates (each of which is internally a vector)

   std::ifstream ist(file.c_str());
   if (ist)
   {
      std::string line_buf;
      while (std::getline(ist, line_buf)) // Read a line in
      {
         std::string sample_name;
         std::stringstream line_stream(line_buf);

         // Read in sample name for the line
         line_stream >> sample_name;

         auto col_it = column_list.begin();
         size_t row = 0;
         while (col_it != column_list.end()) // Go through all required columns
         {
            std::string col_buf; // Work out how many fields to read
            int skip_cols = std::get<0>(*col_it);
            if (col_it != column_list.begin())
            {
               skip_cols -= std::get<0>(*(col_it-1));
            }
            else
            {
               skip_cols -= 1; // Sample name already read
            }

            for (int i = 0; i<skip_cols; ++i)
            {
               line_stream >> col_buf;
            }

            covars_in[row].addCovarValue(sample_name, col_buf);
            ++col_it; // Move to next column
            ++row;
         }
      }
   }
   else
   {
      throw std::runtime_error("Could not open covariate file" + file);
   }

   // Encode categorical columns. Convert strings to doubles
   // All in arma::mat
   arma::mat covar_matrix;
   auto mat_col = covars_in.begin();
   for (auto col_it = column_list.begin(); col_it != column_list.end(); ++col_it)
   {
      // When converting, first sort on sample name so in the same order as
      // pheno vector that will be merged with
      arma::mat new_col;
      if (std::get<1>(*col_it))
      {
         // Categorical covars need dummy encoding
         new_col = encodeDummy(mat_col->sortCovars());
      }
      else
      {
         // Quantitative covars need normalising
         new_col = vecToMat(mat_col->sortCovars());
         normaliseMatCols(new_col);
      }

      // Add the row in
      if (covar_matrix.n_cols > 0)
      {
         covar_matrix = arma::join_rows(covar_matrix, new_col);
      }
      else
      {
         covar_matrix = new_col;
      }
      ++mat_col;
   }

   return covar_matrix;
}

std::vector<std::tuple<int,bool>> parseCovarColumns(const std::string& columns)
{
   std::vector<std::tuple<int,bool>> column_list; // A vector of column numbers, and whether they
                                                  // are categorical and need
                                                  // dummy encoding
   std::string col_buf;
   std::stringstream ss(columns);
   while(std::getline(ss, col_buf, ','))
   {
      if (std::regex_match(col_buf, q_covar_regex))
      {
         auto new_col = std::make_tuple(stoi(std::regex_replace(col_buf, q_covar_regex, "$1")), 0);
         column_list.push_back(new_col);
      }
      else if (std::regex_match(col_buf, covar_regex))
      {
         auto new_col = std::make_tuple(stoi(std::regex_replace(col_buf, covar_regex, "$1")), 1);
         column_list.push_back(new_col);
      }
      else
      {
         std::cerr << "Warning: couldn't parse covariate column " + col_buf << std::endl;
      }
   }

   if (column_list.size() == 0)
   {
      throw std::runtime_error("No covariate columns selected, but file provided");
   }

   // Sort by column number. Use a lambda function
   std::sort(column_list.begin(), column_list.end(),
         [](const std::tuple<int,bool>& a, const std::tuple<int,bool>& b)
         {
            return std::get<0>(a) < std::get<0>(b);
         });

   return column_list;
}

// Converts categorical variates into columns using dummy encoding
arma::mat encodeDummy(const std::vector<std::string>& in_col)
{
   // Find the possible categories, and put them into an associative array (unordered
   // map)
   std::vector<std::string> col_copy = in_col;
   std::sort(col_copy.begin(), col_copy.end());
   auto unique_it = std::unique(col_copy.begin(), col_copy.end());

   std::unordered_map<std::string, int> dummy_map;
   unsigned int dummy_col = 0;
   for (auto cat_it = col_copy.begin(); cat_it != unique_it; ++cat_it)
   {
      dummy_map[*cat_it] = dummy_col++;
   }

   // For each variable, find the category and code this column as 1
   // There will be n-1 categories, as the control is coded as all zeros
   arma::mat out_mat(in_col.size(), dummy_map.size() - 1);
   for (unsigned int i = 0; i < in_col.size(); ++i)
   {
      unsigned int true_col = dummy_map[in_col[i]];
      for (unsigned int j = 0; j < dummy_map.size() - 1; ++j)
      {
         int out_val = 0;
         if (true_col > 0 && true_col - 1 == j)
         {
            out_val = 1;
         }
         out_mat(i, j) = out_val;
      }
   }

   return out_mat;
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

