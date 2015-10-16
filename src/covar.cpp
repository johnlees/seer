/*
 * covar.chpp
 * Helper functions for covar class
 */

#include "covar.hpp"

// Sorts covariates on sample name, returns only values
std::vector<std::string> Covar::sortCovars()
{
   // Sort on sample name
   std::sort(covariates.begin(), covariates.end(),
         [](const std::tuple<std::string, std::string>& lhs, const std::tuple<std::string, std::string>& rhs)
         {
            return std::get<0>(lhs) < std::get<0>(rhs);
         });

   // Return a vector without sample names
   std::vector<std::string> sorted;
   sorted.reserve(covariates.size());
   for (auto it = covariates.begin(); it != covariates.end(); ++it)
   {
      sorted.push_back(std::get<1>(*it));
   }

   return sorted;
}

