/*
 * File: pangwasAssoc.cpp
 *
 * Implements chi^2 and regression association tests for pangwas
 *
 */

#include "pangwas.hpp"

void logisticTest(Kmer& k, const arma::vec& y_train)
{
   // Train classifier
   arma::mat x_train = k.get_x();
   mlpack::regression::LogisticRegression<> fit(x_train.t(), y_train);

   // Extract beta
   arma::vec b_vector = fit.Parameters();
   double b_1 = b_vector(1);
   k.beta(b_1);

   // Extract p-value
   //
   //
   // W = B_1 / SE(B_1) ~ N(0,1)
   //
   // In the special case of a logistic regression, abs can be taken rather
   // than ^2 as responses are 0 or 1
   //
   double W = std::abs(b_1) / pow(varCovarMat(x_train, b_vector)(1,1), 0.5); // null hypothesis b_1 = 0
   double pvalue = normalPval(W);

#ifdef PANGWAS_DEBUG
   std::cerr << "Wald statistic: " << W << "\n";
#endif

   k.p_val(pvalue);
}

//
arma::mat varCovarMat(const arma::mat& x, const arma::mat& b)
{
   // var-covar matrix = inv(I)
   // where I is the Fisher information matrix
   // I = d^2/d(b^2)[log L]
   //
   // see http://czep.net/stat/mlelr.pdf
   arma::mat x_design = join_rows(arma::mat(x.n_rows,1,arma::fill::ones), x);

   // First get logit of x values using parameters from fit, and transform to
   // p(1-p)
   arma::vec y_pred = predictLogitProbs(x_design, b);
   arma::vec y_trans = y_pred % (1 - y_pred);

   // Fill elements of I, which are sums of element by element vector multiples
   arma::mat I(b.n_elem, b.n_elem);
   for (unsigned int i = 0; i<I.n_cols; ++i)
   {
      unsigned int j = 0;
      unsigned int j_max = I.n_rows;

      while (j < j_max)
      {
         I(i,j) = accu(y_trans % x_design.col(i) % x_design.col(j));
         if (i != j)
         {
            I(j,i) = I(i,j); // I is symmetric - only need to calculate upper triangle
         }
         j++;
      }
      j_max--;
   }

   return inv_sympd(I);
}

// returns y = logit(bx)
arma::vec predictLogitProbs(const arma::mat& x, const arma::vec& b)
{
   const arma::vec exponents = x * b;
   const arma::vec y = 1.0 / (1.0 + arma::exp(-exponents));

   return y;
}

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

// Returns p-value for a test statistic that is >0 and standard normally distributed
double normalPval(double testStatistic)
{
   boost::math::normal s;

   double p_val = 2 * (1 - boost::math::cdf(s, testStatistic));
   return p_val;
}
