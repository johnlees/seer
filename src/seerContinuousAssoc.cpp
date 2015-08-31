/*
 * File: seerContinuousAssoc.cpp
 *
 * Implements linear regression association test for seer
 *
 */

#include "seer.hpp"

// Linear fit without covariates
void linearTest(Kmer& k, const arma::vec& y_train)
{
   // Train classifier
   arma::mat x_train = k.get_x();
   doLinear(k, y_train, x_train);
}

// Linear fit with covariates
void linearTest(Kmer& k, const arma::vec& y_train, const arma::mat& mds)
{
   // Train classifier
   arma::mat x_train = arma::join_rows(k.get_x(), mds);
   doLinear(k, y_train, x_train);
}

// Run linear fit
void doLinear(Kmer& k, const arma::vec& y_train, const arma::mat& x_train)
{
   // Singular matrix inversions will throw
   try
   {
      // Set up design matrix, and calculate (X'X)^-1
      arma::mat x_design = join_rows(arma::mat(x_train.n_rows,1,arma::fill::ones), x_train);
      arma::mat inv_xxt = inv_sympd(x_design*x_design.t());

      // Regress beta: b = (X'X)^-1*Xy
      arma::vec b = inv_xxt * x_design * y_train;
      k.beta(b(1));

      // Extract p-value
      //
      // W = B_1 / SE(B_1) ~ N(0,1)
      //
      // SE(B_1) = MSE * (X'X)^-1
      // MSE = sum(Y_i-Y'_i)^2 / n-2
      //
      double MSE = accu(square(y_train - predictLinearProbs(x_design, b))) / (x_train.n_rows - 2);
      double W = std::abs(k.beta()) / (inv_xxt(1,1)*MSE); // null hypothesis b_1 = 0
      k.p_val(normalPval(W));

#ifdef SEER_DEBUG
      std::cerr << "Wald statistic: " << W << "\n";
      std::cerr << "p-value: " << k.p_val() << "\n";
#endif
   }
   // Method will throw if a singular matrix is inverted
   catch (std::exception& e)
   {
      std::cerr << k.sequence() << "\n"
                << "kmer convergence error: "
                << e.what() << "\n";

      k.p_val(0);
      k.beta(0);
   }
}

// returns y' = bx
arma::vec predictLinearProbs(const arma::mat& x, const arma::vec& b)
{
   const arma::vec y = x * b;

   return y;
}

