/*
 * File: pangwasAssoc.cpp
 *
 * Implements chi^2 and regression association tests for pangwas
 *
 */

#include "pancommon.hpp"

// Basic chi^2 test, using contingency table
double chiTest(arma::mat& table)
{
   double chisq = 0;
   int N = accu(table);

   if (N == 0)
   {
      throw std::logic_error("Empty table for chisq test\n");
   }

   // Without Yates' continuity correction
   chisq = N * pow(det(table), 2);
   for (int i = 0; i < 2; ++i)
   {
      chisq /= accu(table.row(i)) * accu(table.col(i));
   }

   // For df = 1, as here, chi^2 == N(0,1)^2 (standard normal dist.)
   double p_value = normalPval(pow(chisq, 0.5));
#ifdef PANGWAS_DEBUG
   std::cerr << "chisq:" << chisq << "\n";
   std::cerr << "chisq p: " << p_value << "\n";
#endif
   return p_value;
}

// Welch two sample t-test, for continuous phenotypes
double welchTwoSamplet(const arma::vec& x, const arma::vec& y)
{
   // Subset into present and absent groups
   arma::vec group1 = y.elem(find(x==0));
   arma::vec group2 = y.elem(find(x==1));

   // Calculate group means and variances
   double x1 = mean(group1);
   double x2 = mean(group2);
   double v1 = var(group1);
   double v2 = var(group2);

   // t and degrees freedom for test
   double t = (x1 - x2)*pow((v1/group1.n_rows + v2/group2.n_rows), -0.5);
   double df = pow((v1/group1.n_rows + v2/group2.n_rows), 2) / (pow(v1/group1.n_rows,2)/(group1.n_rows-1) + pow(v2/group2.n_rows,2)/(group2.n_rows-1));

   // Calculate p-value from t distribution
   boost::math::students_t t_dist(df);
   double p_val = 2 * (1 - boost::math::cdf(t_dist, t));
#ifdef PANGWAS_DEBUG
   std::cerr << "welch t:" << t << "df:" << df << "\n";
   std::cerr << "welch p-val:" << p_val << "\n";
#endif

   return p_val;
}

// Returns p-value for a test statistic that is >0 and standard normally distributed
double normalPval(double testStatistic)
{
   boost::math::normal s;

   double p_val = 2 * (1 - boost::math::cdf(s, testStatistic));
   return p_val;
}
