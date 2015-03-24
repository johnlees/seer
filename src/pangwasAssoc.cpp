/*
 * File: pangwasAssoc.cpp
 *
 * Implements chi^2 and regression association tests for pangwas
 *
 */

#include "pangwas.h"

void logisticTest(Kmer& k, arma::vec y_train)
{
   // Train classifier
   arma::vec x_train = k.get_x();
   mlpack::regression::LogisticRegression<> fit(x_train, y_train);

   // Extract beta
   double b_1 = fit.Parameters()[1];
   k.beta(b_1);

   // Extract p-value
   arma::vec y_predicted;
   fit.Predict(x_train, y_predicted);

   //
   // SE^2 = sum[(Y_i - Y'_i)^2]
   //        -------------------
   //               N - 2
   //
   // W = B_1 / SE(B_1) ~ N(0,1)
   //
   // In the special case of a logistic regression, abs can be taken rather
   // than ^2 as responses are 0 or 1
   //
   double b1_se = pow(accu(abs(y_train - y_predicted)) / (y_train.n_rows - 2), 0.5);
   double W = b_1/b1_se; // b_0 = 0

   double pvalue = normalPval(W);
   k.p_val(pvalue);
}

double chiTest(arma::Mat<int>& table)
{
   // Use floats, speed over accuracy?
   double chisq = 0;
   int N = accu(table);

   if (N == 0)
   {
      throw std::logic_error("Empty table for chisq test\n");
   }

   for (int i = 0; i < 2; ++i)
   {
      for (int j = 0; j < 2; ++j)
      {
         double fe = (accu(table.row(i)) * accu(table.col(j))) / N;
         if (fe == 0)
         {
            throw std::runtime_error("Table element " + std::to_string(i) + "'" + std::to_string(j) + " fe = 0\n");
         }

         chisq += pow(table(i,j) - fe, 2);
      }
   }

   // For df = 1, as here, chi^2 == N(0,1)^2
   double p_value = normalPval(pow(chisq, 0.5));
#ifdef PANGWAS_DEBUG
   std::cerr << "chisq result: " << p_value << "\n";
#endif
   return p_value;
}

// Returns p-value for a test statistic that is >0 and standard normally distributed
double normalPval(double testStatistic)
{
   boost::math::normal s;

   double p_val = 2 * (1 - boost::math::cdf(s, testStatistic));
   return p_val;
}
