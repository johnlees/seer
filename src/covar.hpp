/*
 * covar.hpp
 * Header file for covar class
 */

#include <string>
#include <vector>
#include <tuple>
#include <algorithm>

// Samples have a name and phenotype
// Covariates will be added later
class Covar
{
   public:
      // Initialisation
      Covar() {}

      // Add values, and read back
      std::vector<std::string> sortCovars();
      void addCovarValue (const std::string& sample, const std::string& value)
         {covariates.push_back(std::make_tuple(sample, value));}

   private:
      std::vector<std::tuple<std::string, std::string>> covariates;
};

