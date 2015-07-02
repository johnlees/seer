/*
 * File: pangwasBinaryAssoc.cpp
 *
 * Implements logistic regression association tests for pangwas
 *
 */

#include "pangwas.hpp"

// Logistic fit without covariates
void logisticTest(Kmer& k, const arma::vec& y_train, const unsigned int nr)
{
   // Train classifier
   arma::mat x_train = k.get_x();
   doLogit(k, y_train, x_train, nr);
}

// Logistic fit with covariates
void logisticTest(Kmer& k, const arma::vec& y_train, const unsigned int nr, const arma::mat& mds)
{
   // Train classifier
   arma::mat x_train = arma::join_rows(k.get_x(), mds);
   doLogit(k, y_train, x_train, nr);
}

void doLogit(Kmer& k, const arma::vec& y_train, const arma::mat& x_train, const unsigned int nr)
{
   regression fit;
   try
   {
      if (nr != 1)
      {
         fit = logisticPval(y_train, x_train);
      }
      else
      {
         fit = newtonRaphson(y_train, x_train);
      }
   }
   // Methods will throw if a singular matrix is inverted
   catch (std::exception& e)
   {
      std::cerr << k.sequence() << "\n"
                << "kmer convergence error: "
                << e.what() << "\n";

      fit.p_val = 0;
      fit.beta = 0;
   }

   k.p_val(fit.p_val);
   k.beta(fit.beta);
}

regression newtonRaphson(const arma::vec& y_train, const arma::mat& x_train)
{
   regression parameters;

   // Keep iterations to track convergence
   // Also useful to keep second derivative, for calculating p-value
   std::vector<arma::vec> parameter_iterations;
   arma::mat var_covar_mat;

   // Get starting point from a linear regression, which is fast
   // and will reduce number of n-r iterations
   // Set up design matrix, and calculate (X'X)^-1
   arma::mat x_design = join_rows(arma::mat(x_train.n_rows,1,arma::fill::ones), x_train);
   try
   {
      arma::mat inv_xxt = inv_sympd(x_design*x_design.t());

      // Regress beta: b = (X'X)^-1*Xy
      parameter_iterations.push_back(inv_xxt * x_design * y_train);
   }
   // Method will throw if a singular matrix is inverted
   catch (std::exception& e)
   {
      // Alternatively just use:
      parameter_iterations.push_back(arma::ones(x_train.n_cols + 1));
   }

   for (unsigned int i = 0; i < max_nr_iterations; ++i)
   {
      arma::vec b0 = parameter_iterations.back();
      arma::vec y_pred = predictLogitProbs(x_design, b0);

      var_covar_mat = inv_sympd(x_design.t() * diagmat(y_pred % (arma::ones(y_pred.n_rows) - y_pred)) * x_design);
      arma::vec b1 = b0 + var_covar_mat * x_design.t() * (y_train - y_pred);
      parameter_iterations.push_back(b1);

      if (std::abs(b1(1) - b0(1)) < convergence_limit)
      {
         break;
      }
   }

#ifdef PANGWAS_DEBUG
   std::cerr << "Number of iterations: " << parameter_iterations.size() << "\n";
#endif

   parameters.beta = parameter_iterations.back()(1);

   double W = std::abs(parameters.beta) / pow(var_covar_mat(1,1), 0.5);
   parameters.p_val = normalPval(W);

 #ifdef PANGWAS_DEBUG
   std::cerr << "Wald statistic: " << W << "\n";
   std::cerr << "p-value: " << parameters.p_val << "\n";
#endif

   return parameters;
}

// Fits logistic regression, and returns beta and pvalue
regression logisticPval(const arma::vec& y_train, const arma::mat& x_train)
{
   regression parameters;

   arma::mat x_design = join_rows(arma::mat(x_train.n_rows,1,arma::fill::ones), x_train);
   column_vector starting_point(x_design.n_cols);
/*   try
   {
      arma::mat inv_xxt = inv_sympd(x_design*x_design.t());

      // Regress beta: b = (X'X)^-1*Xy
      starting_point = arma_to_dlib(inv_xxt * x_design * y_train);
   }
   // Method will throw if a singular matrix is inverted
   catch (std::exception& e)
   {
      // Alternatively just use:
      for (size_t i = 0; i < x_design.n_cols; ++i)
      {
         starting_point(i) = 1;
      }
   }
*/
   for (size_t i = 0; i < x_design.n_cols; ++i)
   {
      starting_point(i) = 1;
   }

   // Use BFGS optimiser in dlib to maximise likelihood function by chaging the
   // b vector, which will end in starting_point
   dlib::find_max(dlib::bfgs_search_strategy(),
                  dlib::objective_delta_stop_strategy(convergence_limit),
                  LogitLikelihood(x_train, y_train), LogitLikelihoodGradient(x_train, y_train),
                  starting_point, -1);

   // Extract beta
   arma::vec b_vector = dlib_to_arma(starting_point);
   double b_1 = b_vector(1);
   parameters.beta = b_1;

   // Extract p-value
   //
   //
   // W = B_1 / SE(B_1) ~ N(0,1)
   //
   // In the special case of a logistic regression, abs can be taken rather
   // than ^2 as responses are 0 or 1
   //
   double W = std::abs(b_1) / pow(varCovarMat(x_design, b_vector)(1,1), 0.5); // null hypothesis b_1 = 0
   parameters.p_val = normalPval(W);

#ifdef PANGWAS_DEBUG
   std::cerr << "Wald statistic: " << W << "\n";
   std::cerr << "p-value: " << parameters.p_val << "\n";
#endif

   return parameters;
}

// Returns var-covar matrix for logistic function
// WARNING: This contains an inversion, which will throw on a singular matrix
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

   return inv_sympd(I);
}

// returns y = logit(bx)
arma::vec predictLogitProbs(const arma::mat& x, const arma::vec& b)
{
   const arma::vec exponents = x * b;
   const arma::vec y = 1.0 / (1.0 + arma::exp(-exponents));

   return y;
}
