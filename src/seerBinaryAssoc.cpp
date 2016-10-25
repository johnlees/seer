/*
 * File: seerBinaryAssoc.cpp
 *
 * Implements logistic regression association tests for seer
 *
 */

#include "seer.hpp"

// Logistic fit without covariates
void logisticTest(Kmer& k, const arma::vec& y_train, const double null_ll)
{
   // Train classifier
   arma::mat x_train = arma::join_rows(arma::mat(y_train.n_rows,1,arma::fill::ones), k.get_x());
   doLogit(k, y_train, x_train);

   // Likelihood ratio test
   k.lrt_p_val(likelihoodRatioTest(k, null_ll));
}

// Logistic fit with covariates
void logisticTest(Kmer& k, const arma::vec& y_train, const double null_ll, const arma::mat& mds)
{
   // Train classifier
   arma::mat x_train = arma::join_rows(arma::mat(y_train.n_rows,1,arma::fill::ones), k.get_x());
   x_train = arma::join_rows(x_train, mds);
   doLogit(k, y_train, x_train);

   // Likelihood ratio test
   k.lrt_p_val(likelihoodRatioTest(k, null_ll));
}

// This uses BFGS optimisation by default. Invokes NR or Firth on error
void doLogit(Kmer& k, const arma::vec& y_train, const arma::mat& x_design)
{
   column_vector starting_point(x_design.n_cols);

   if (k.firth())
   {
      newtonRaphson(k, y_train, x_design, 1);
   }
   else
   {
      starting_point(0) = log(mean(y_train)/(1 - mean(y_train)));
      for (size_t i = 1; i < x_design.n_cols; ++i)
      {
         starting_point(i) = bfgs_start_beta;
      }

      try
      {
         // Use BFGS optimiser in dlib to maximise likelihood function by chaging the
         // b vector, which will end in starting_point
         LogitLikelihood likelihood_fit(x_design, y_train); // store this, as it is used for computing the LRT

         dlib::find_max(dlib::bfgs_search_strategy(),
                     dlib::objective_delta_stop_strategy(convergence_limit),
                     likelihood_fit, LogitLikelihoodGradient(x_design, y_train),
                     starting_point, -1);

         // Extract beta and likelihood
         arma::vec b_vector = dlib_to_arma(starting_point);
         k.beta(b_vector(1));

         k.log_likelihood(likelihood_fit(starting_point));

         // Extract p-value
         //
         //
         // W = B_1 / SE(B_1) ~ N(0,1)
         //
         // In the special case of a logistic regression, abs can be taken rather
         // than ^2 as responses are 0 or 1
         //
         arma::mat var_covar_mat = varCovarMat(x_design, b_vector);
         double se = pow(var_covar_mat(1,1), 0.5);

         // Zeros will result in bad regression with large SE - firth regression helps
         if (se > se_limit)
         {
            throw std::runtime_error("se>limit");
         }
         else
         {
            k.standard_error(se);

            double W = std::abs(b_vector(1)) / se; // null hypothesis b_1 = 0
            k.p_val(normalPval(W));

#ifdef SEER_DEBUG
            std::cerr << "Wald statistic: " << W << "\n";
            std::cerr << "p-value: " << k.p_val() << "\n";
#endif
            // Add in covariate p-values
            for (unsigned int i = 2; i < var_covar_mat.n_rows; ++i)
            {
               se = pow(var_covar_mat(i,i), 0.5);
               W = std::abs(b_vector(i)) / se;

               k.add_covar_p(normalPval(W));
            }
         }
      }
      // Sometimes won't converge, use N-R instead
      catch (std::exception& e)
      {
#ifdef SEER_DEBUG
         std::cerr << "Caught error " << e.what() << std::endl;
#endif

         // SE is greater than specified limit - run Firth regression
         if (strcmp(e.what(), "se>limit") == 0)
         {
            k.add_comment("large-se");
            newtonRaphson(k, y_train, x_design, 1);
         }
         // BFGS optimiser did not converge - use NR iterations w/o Firth first
         // Could also be matrix inversion failing
         else
         {
            k.add_comment("bfgs-fail");
            newtonRaphson(k, y_train, x_design);
         }
      }
   }
}

void newtonRaphson(Kmer& k, const arma::vec& y_train, const arma::mat& x_design, const bool firth)
{
   // Keep iterations to track convergence
   // Also useful to keep second derivative, for calculating p-value
   std::vector<arma::vec> parameter_iterations;
   arma::mat var_covar_mat;
   int failed = 0;

   // Could get starting point from a linear regression, which is fast
   // and will reduce number of n-r iterations
   // Set up design matrix, and calculate (X'X)^-1
   // Seems more reliable to go for b = 0, plus a non-zero intercept
   // See: doi:10.1016/S0169-2607(02)00088-3
   arma::vec starting_point = arma::zeros(x_design.n_cols);
   starting_point(0) = log(mean(y_train)/(1 - mean(y_train)));

   parameter_iterations.push_back(starting_point);

   for (unsigned int i = 0; i < max_nr_iterations; ++i)
   {
      arma::vec b0 = parameter_iterations.back();
      arma::vec y_pred = predictLogitProbs(x_design, b0);

      arma::mat U(x_design.n_cols, 1);
      arma::mat W = repmat(y_pred % (arma::ones(y_pred.n_rows) - y_pred), 1, x_design.n_cols);

      // Perform inversion, which may fail
      var_covar_mat = inv_covar(x_design.t() * (W % x_design));
      if (var_covar_mat.n_cols == 0 || var_covar_mat.n_rows == 0)
      {
         k.add_comment("inv-fail");
         k.p_val(0);
         std::cerr << "Inversion at input line " << k.line_number() << " failed" << std::endl;
         failed = 1;
         break;
      }

      if (firth)
      {
         // Firth logistic regression
         // See: DOI: 10.1002/sim.1047
         // Hat matrix
         // Note: W is diagonal so X.t() * W * X is still sympd
         arma::mat H = (sqrt(W) % x_design) * var_covar_mat * (x_design.t() % sqrt(W).t());

         arma::vec correction(y_train.n_rows);
         correction.fill(0.5);

         // Penalised score
         U = x_design.t() * (y_train - y_pred + diagvec(H) % (correction - y_pred));
      }
      else
      {
         U = x_design.t() * (y_train - y_pred);
      }

      arma::vec b1 = b0 + var_covar_mat * U;
      parameter_iterations.push_back(b1);

      if (std::abs(b1(1) - b0(1)) < convergence_limit)
      {
         break;
      }
   }

#ifdef SEER_DEBUG
   std::cerr << "Number of iterations: " << parameter_iterations.size() - 1 << "\n";
#endif
   // If convergence not reached, try Firth logistic regression
   if (parameter_iterations.size() == max_nr_iterations)
   {
      if (!firth)
      {
         k.add_comment("nr-fail");
         newtonRaphson(k, y_train, x_design, 1);
      }
      else
      {
         k.add_comment("firth-fail");
      }
   }
   else if (!failed)
   {
      // Add beta and log-likelihood
      column_vector converged_beta = arma_to_dlib(parameter_iterations.back());

      LogitLikelihood likelihood_fit(x_design, y_train);
      if (k.firth())
      {
         k.log_likelihood(likelihood_fit(converged_beta) + 0.5*log(det(inv_covar(var_covar_mat))));
      }
      else
      {
         k.log_likelihood(likelihood_fit(converged_beta));
      }

      k.beta(converged_beta(1));

      double se = pow(var_covar_mat(1,1), 0.5);
      k.standard_error(se);

      // Deal with large SEs
      if (se > se_limit)
      {
         if (!firth)
         {
            newtonRaphson(k, y_train, x_design, 1);
         }
         else
         {
            k.add_comment("large-se");
         }
      }

      double W = std::abs(k.beta()) / se;
      k.p_val(normalPval(W));

#ifdef SEER_DEBUG
      std::cerr << "Wald statistic: " << W << "\n";
      std::cerr << "p-value: " << k.p_val() << "\n";
#endif
      // Add in covariate p-values
      for (unsigned int i = 2; i< var_covar_mat.n_rows; ++i)
      {
         se = pow(var_covar_mat(i,i), 0.5);
         W = std::abs(parameter_iterations.back()(i)) / se;

         k.add_covar_p(normalPval(W));
      }
   }
}

// Returns var-covar matrix for logistic function
arma::mat varCovarMat(const arma::mat& x, const arma::mat& b)
{
   // var-covar matrix = inv(I)
   // where I is the Fisher information matrix
   // I = d^2/d(b^2)[log L]
   //
   // see http://czep.net/stat/mlelr.pdf

   // First get logit of x values using parameters from fit, and transform to
   // p(1-p)
   arma::vec y_pred = predictLogitProbs(x, b);
   arma::vec y_trans = y_pred % (1 - y_pred);

   // Fill elements of I, which are sums of element by element vector multiples
   arma::mat I(b.n_elem, b.n_elem);
   unsigned int j_max = I.n_rows;
   for (unsigned int i = 0; i<I.n_cols; ++i)
   {
      for (unsigned int j = i; j < j_max; j++)
      {
         I(i,j) = accu(y_trans % x.col(i) % x.col(j));
         if (i != j)
         {
            I(j,i) = I(i,j); // I is symmetric - only need to calculate upper triangle
         }
      }
   }

   return inv_covar(I);
}

// returns y = logit(bx)
arma::vec predictLogitProbs(const arma::mat& x, const arma::vec& b)
{
   const arma::vec exponents = x * b;
   const arma::vec y = 1.0 / (1.0 + arma::exp(-exponents));

   return y;
}
