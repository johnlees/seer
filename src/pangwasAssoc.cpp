/*
 * File: pangwasAssoc.cpp
 *
 * Implements chi^2 and regression association tests for pangwas
 *
 */

#include "pangwas.h"

void logisticTest(Kmer& k, arma::vec y, const double p_cutoff)
{
   // Train classifier
   mlpack::regression::LogisticRegression<> fit(k.get_x(), y);

   // TODO
   // Extract pvalue
   double pvalue = 1;
   double beta = fit.Parameters()[0]; //TODO correct? Transpose needed?

   k.p_val(pvalue);
   k.beta(beta);
}

float chiTest(arma::Mat<int>& table)
{
   // Use floats, speed over accuracy?
   float chisq = 0;
   int N = accu(table);

   if (N == 0)
   {
      throw std::logic_error("Empty table for chisq test\n");
   }

   for (int i = 0; i < 2; ++i)
   {
      for (int j = 0; j < 2; ++j)
      {
         float fe = (accu(table.row(i)) * accu(table.col(j))) / N;
         if (fe == 0)
         {
            throw std::runtime_error("Table element " + std::to_string(i) + "'" + std::to_string(j) + " fe = 0\n");
         }

         chisq += pow(table(i,j) - fe, 2);
      }
   }

   // For df = 1, as here, chi^2 == N(0,1)^2
   boost::math::normal s;
   float p_val = 2 * (1 - boost::math::cdf(s, pow(chisq,0.5)));

   return p_val;
}
