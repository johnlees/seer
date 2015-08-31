/*
 * logitFunction.hpp
 * Header file for logitFunction class
*/

#include "seercommon.hpp"

class LogitFunction
{
   public:
      // Initialisation
      LogitFunction(const arma::mat& _predictors, const arma::vec& _responses, const double _lambda = 0);

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

class LogitLikelihood : public LogitFunction
{
   public:
      LogitLikelihood(const arma::mat& _predictors, const arma::vec& _responses, const double _lambda = 0)
         : LogitFunction(_predictors, _responses, _lambda)
      {
      }

      double operator() (const column_vector& parameters_in) const;
};

class LogitLikelihoodGradient : public LogitFunction
{
   public:
      LogitLikelihoodGradient(const arma::mat& _predictors, const arma::vec& _responses, const double _lambda = 0)
         : LogitFunction(_predictors, _responses, _lambda)
      {
      }

      column_vector operator() (const column_vector& parameters_in) const;
};
