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
   // Set up design matrix
   arma::mat x_design = join_rows(arma::mat(x_train.n_rows,1,arma::fill::ones), x_train);
   arma::vec b;

   try
   {
      // A starting point for the optimiser
      column_vector starting_point(x_design.n_cols);
      starting_point(0) = mean(y_train);
      for (size_t i = 1; i < x_design.n_cols; ++i)
      {
         starting_point(i) = bfgs_start_beta;
      }

      // Use BFGS optimiser in dlib to minimise OLS difference by chaging the
      // b vector, which will end in starting_point
      dlib::find_min(dlib::bfgs_search_strategy(),
                  dlib::objective_delta_stop_strategy(convergence_limit),
                  LinearLikelihood(x_design, y_train), LinearLikelihoodGradient(x_design, y_train),
                  starting_point, -1);

      // Extract beta
      b = dlib_to_arma(starting_point);

   }
   catch (std::exception& e)
   {
#ifdef SEER_DEBUG
      std::cerr << "bfgs failed with " << e.what() << std::endl;
#endif
      k.add_comment("bfgs-fail");

      // Calculate (X'X)^-1. Use QR decomposition to solve.
      // Might be slower than Cholesky, but numerically is more stable
      arma::mat Q, R;

      // Regress beta: b = (X'X)^-1*Xy
      // i.e. solve Rb = Q.t() * y where X=QR
      arma::qr(Q, R, x_design);
      arma::solve(b, R, Q.t() * y_train);
   }

   k.beta(b(1));

   // Extract p-value
   //
   // W = B_1 / SE(B_1) ~ N(0,1)
   //
   // SE(B_1) = MSE * (X'X)^-1
   // MSE = sum(Y_i-Y'_i)^2 / n-2
   //
   double MSE = accu(square(y_train - predictLinearProbs(x_design, b))) / (x_train.n_rows - 2);
   arma::mat var_covar_mat = inv_covar(x_design.t()*x_design);

   double se = pow((var_covar_mat(1,1) * MSE), 0.5);
   k.standard_error(se);

   double W = std::abs(k.beta()) / (se); // null hypothesis b_1 = 0
   k.p_val(normalPval(W));

#ifdef SEER_DEBUG
   std::cerr << "Wald statistic: " << W << "\n";
   std::cerr << "p-value: " << k.p_val() << "\n";
#endif

   // Add in covariate p-values
   for (unsigned int i = 2; i < var_covar_mat.n_rows; ++i)
   {
      se = pow((var_covar_mat(i,i) * MSE), 0.5);
      W = std::abs(b(i)) / (se);

      k.add_covar_p(normalPval(W));
   }
}

// returns y' = bx
arma::vec predictLinearProbs(const arma::mat& x, const arma::vec& b)
{
   const arma::vec y = x * b;

   return y;
}

