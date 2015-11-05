/*
 * File: linearFunction.cpp
 *
 * Likelihood, gradient and initialisation of linear function
 * Functor, based on mlpack: www.mlpack.org
 */

#include "linkFunction.hpp"

/**
 * Evaluate the logistic regression objective function given the estimated
 * parameters.
 */
double LinearLikelihood::operator()(const column_vector& parameters_in)
   const
{
   // Convert from dlib column matrix to armadillo column matrix
   arma::vec parameters = dlib_to_arma(parameters_in);

   // L(b) = 1/2*||y-Xb||^2
   double result = 0.5 * accu(square(responses - predictors * parameters));

   return result;
}

// Evaluate the gradient of the logistic regression objective function.
column_vector LinearLikelihoodGradient::operator()(const column_vector& parameters_in)
   const
{
   // Convert from dlib column matrix to armadillo column matrix
   arma::vec parameters = dlib_to_arma(parameters_in);
   arma::vec gradient(parameters.n_elem);

   // dL(b)/db = X.t()(Xb - y)
   gradient = predictors.t() * (predictors * parameters - responses);

   return arma_to_dlib(gradient);
}

