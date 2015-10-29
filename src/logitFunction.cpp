/*
 * File: logitFunction.cpp
 *
 * Likelihood, gradient and initialisation of logit function
 * Functor, based on mlpack: www.mlpack.org
 */

#include "logitFunction.hpp"

LogitFunction::LogitFunction(const arma::mat& _predictors, const arma::vec& _responses, const double _lambda)
   :predictors(_predictors), responses(_responses), lambda(_lambda)
{
}

/**
 * Evaluate the logistic regression objective function given the estimated
 * parameters.
 */
double LogitLikelihood::operator()(const column_vector& parameters_in)
    const
{
   // Convert from dlib column matrix to armadillo column matrix
   arma::vec parameters = dlib_to_arma(parameters_in);

   // The objective function is the log-likelihood function (w is the parameters
   // vector for the model; y is the responses; x is the predictors; sig() is the
   // sigmoid function):
   //   f(w) = sum(y log(sig(w'x)) + (1 - y) log(sig(1 - w'x))).
   // We want to minimize this function.  L2-regularization is just lambda
   // multiplied by the squared l2-norm of the parameters then divided by two.

   // For the regularization, we ignore the first term, which is the intercept
   // term.
   const double regularization = 0.5 * lambda *
       arma::dot(parameters.col(0).subvec(1, parameters.n_elem - 1),
                 parameters.col(0).subvec(1, parameters.n_elem - 1));

   // Calculate vectors of sigmoids.  The intercept term is parameters(0, 0) and
   // does not need to be multiplied by any of the predictors.
   const arma::vec exponents = parameters(0) + predictors *
       parameters.col(0).subvec(1, parameters.n_elem - 1);
   const arma::vec sigmoid = 1.0 / (1.0 + arma::exp(-exponents));

   // Assemble full objective function.  Often the objective function and the
   // regularization as given are divided by the number of features, but this
   // doesn't actually affect the optimization result, so we'll just ignore those
   // terms for computational efficiency.
   double result = 0.0;
   for (size_t i = 0; i < responses.n_elem; ++i)
   {
     if (responses[i] == 1)
       result += log(sigmoid[i]);
     else
       result += log(1.0 - sigmoid[i]);
   }

   return result - regularization;
}

// Evaluate the gradient of the logistic regression objective function.
column_vector LogitLikelihoodGradient::operator()(const column_vector& parameters_in)
   const
{
   // Convert from dlib column matrix to armadillo column matrix
   arma::vec parameters = dlib_to_arma(parameters_in);
   arma::vec gradient(parameters.n_elem);

   // Regularization term.
   arma::mat regularization;
   regularization = lambda * parameters.col(0).subvec(1, parameters.n_elem - 1);

   const arma::vec sigmoids = 1 / (1 + arma::exp(-parameters(0, 0)
       - predictors * parameters.col(0).subvec(1, parameters.n_elem - 1)));

   gradient[0] = arma::accu(responses - sigmoids);
   gradient.col(0).subvec(1, parameters.n_elem - 1) = predictors.t() * (responses -
       sigmoids) - regularization;

   return arma_to_dlib(gradient);
}

