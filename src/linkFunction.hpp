/*
 * linkFunction.hpp
 * Header file for logitFunction and linearFunction class
*/

#include "seercommon.hpp"

class LinkFunction
{
   public:
      // Initialisation
      LinkFunction(const arma::mat& _predictors, const arma::vec& _responses, const double _lambda = 0)
         : predictors(_predictors), responses(_responses), lambda(_lambda)
      {
      }

      // Likelihood and first derivative
      double likelihood(const column_vector& b)
         const;
      column_vector gradient(const column_vector& b)
         const;

   protected:
      arma::mat predictors;
      arma::vec responses;

      double lambda;
};

class LogitLikelihood : public LinkFunction
{
   public:
      LogitLikelihood(const arma::mat& _predictors, const arma::vec& _responses, const double _lambda = 0)
         : LinkFunction(_predictors, _responses, _lambda)
      {
      }

      double operator() (const column_vector& parameters_in) const;
};

class LogitLikelihoodGradient : public LinkFunction
{
   public:
      LogitLikelihoodGradient(const arma::mat& _predictors, const arma::vec& _responses, const double _lambda = 0)
         : LinkFunction(_predictors, _responses, _lambda)
      {
      }

      column_vector operator() (const column_vector& parameters_in) const;
};

class LinearLikelihood : public LinkFunction
{
   public:
      LinearLikelihood(const arma::mat& _predictors, const arma::vec& _responses, const double _lambda = 0)
         : LinkFunction(_predictors, _responses, _lambda)
      {
      }

      double operator() (const column_vector& parameters_in) const;
};

class LinearLikelihoodGradient : public LinkFunction
{
   public:
      LinearLikelihoodGradient(const arma::mat& _predictors, const arma::vec& _responses, const double _lambda = 0)
         : LinkFunction(_predictors, _responses, _lambda)
      {
      }

      column_vector operator() (const column_vector& parameters_in) const;
};
