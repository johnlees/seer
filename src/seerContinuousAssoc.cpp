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
   // Set up design matrix, and calculate (X'X)^-1
   // Use QR decomposition to solve. Might be slower than Cholesky, but
   // numerically is more stable
   arma::mat Q, R;
   arma::vec b;
   arma::mat x_design = join_rows(arma::mat(x_train.n_rows,1,arma::fill::ones), x_train);

   // Regress beta: b = (X'X)^-1*Xy
   // i.e. solve Rb = Q.t() * y where X=QR
   arma::qr(Q, R, x_design);
   arma::solve(b, R, Q.t() * y_train);
   k.beta(b(1));

   // Extract p-value
   //
   // W = B_1 / SE(B_1) ~ N(0,1)
   //
   // SE(B_1) = MSE * (X'X)^-1
   // MSE = sum(Y_i-Y'_i)^2 / n-2
   //
   double MSE = accu(square(y_train - predictLinearProbs(x_design, b))) / (x_train.n_rows - 2);
   double se = pow((inv_covar(x_design.t()*x_design)(1,1) * MSE), 0.5);
   k.standard_error(se);

   double W = std::abs(k.beta()) / (se); // null hypothesis b_1 = 0
   k.p_val(normalPval(W));

#ifdef SEER_DEBUG
      std::cerr << "Wald statistic: " << W << "\n";
      std::cerr << "p-value: " << k.p_val() << "\n";
#endif
}

// returns y' = bx
arma::vec predictLinearProbs(const arma::mat& x, const arma::vec& b)
{
   const arma::vec y = x * b;

   return y;
}

