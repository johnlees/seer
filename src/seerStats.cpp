/*
 * File: seerChiFilter
 *
 * Implements chi^2 and welch two sample t test for seer and kmds filtering
 *
 */

#include "seer.hpp"

const double normalArea = pow(2*M_PI, -0.5);

// Basic chi^2 test, using contingency table
double chiTest(Kmer& k, const arma::vec& y)
{
   arma::mat x = k.get_x();

   double chisq = 0;

   // Contigency table
   //         unaffected affected
   // present a          b
   // absent  c          d
   //
   // Use doubles for compatibility with det function in arma::mat
   double a = 0, b = 0, c = 0, d = 0;

   arma::vec::const_iterator j = y.begin();
   for (arma::vec::const_iterator i = x.begin(); i!=x.end(); ++i)
   {
      if (*j == 0) {
         if (*i == 0){
            c++;
         } else {
            a++;
         }
      } else {
         if (*i == 0){
            d++;
         } else {
            b++;
         }
      }
      j++;
   }

   arma::mat::fixed<2, 2> table = {a, b, c, d};
#ifdef SEER_DEBUG
   arma::Mat<int>::fixed<2, 2> tab_out = {int (a), int (b), int (c), int(d)};
   std::cerr << tab_out << "\n";
#endif

   int N = accu(table);

   if (N == 0)
   {
      throw std::logic_error("Empty table for chisq test\n");
   }

   // Treat as invalid if any entry is 0 or 1, or if more than one entry < 5
   // Mark as needing to use Firth regression
   int low_obs = 0;
   for (auto obs = table.begin(); obs != table.end(); ++obs)
   {
      if (*obs <= 1 || (*obs <= 5 && ++low_obs > 2))
      {
         k.add_comment("bad-chisq");
         k.firth(1);
         break;
      }
   }

   // Without Yates' continuity correction
   chisq = N * pow(det(table), 2);
   for (int i = 0; i < 2; ++i)
   {
      chisq /= accu(table.row(i)) * accu(table.col(i));
   }

   // For df = 1, as here, chi^2 == N(0,1)^2 (standard normal dist.)
   double p_value = normalPval(pow(chisq, 0.5));
#ifdef SEER_DEBUG
   std::cerr << "chisq:" << chisq << "\n";
   std::cerr << "chisq p: " << p_value << "\n";
#endif
   return p_value;
}

// Welch two sample t-test, for continuous phenotypes
double welchTwoSamplet(const Kmer& k, const arma::vec& y)
{
   arma::mat x = k.get_x();

   // Subset into present and absent groups
   arma::vec group1 = y.elem(find(x==0));
   arma::vec group2 = y.elem(find(x==1));

   // Calculate group means and variances
   double p_val = 0;
   if (group1.n_elem != 0 && group2.n_elem != 0)
   {
      double x1 = mean(group1);
      double x2 = mean(group2);
      double v1 = var(group1);
      double v2 = var(group2);

      // t and degrees freedom for test
      double t = (x1 - x2)*pow((v1/group1.n_rows + v2/group2.n_rows), -0.5);
      double df = pow((v1/group1.n_rows + v2/group2.n_rows), 2) / (pow(v1/group1.n_rows,2)/(group1.n_rows-1) + pow(v2/group2.n_rows,2)/(group2.n_rows-1));

      // Calculate p-value from t distribution
      boost::math::students_t t_dist(df);
      p_val = 2 * (1 - boost::math::cdf(t_dist, t));
#ifdef SEER_DEBUG
      std::cerr << "welch t:" << t << "df:" << df << "\n";
      std::cerr << "welch p-val:" << p_val << "\n";
#endif
   }

   return p_val;
}

// Fit null models for null log-likelihoods
double nullLogLikelihood(const arma::mat& x, const arma::vec& y, const int continuous)
{
   double null_ll = 0;
   Kmer null_kmer;
   if (x.n_cols > 1)
   {
      if (continuous)
      {
         doLinear(null_kmer, y, x);
      }
      else
      {
         doLogit(null_kmer, y, x);
      }
      null_ll = null_kmer.log_likelihood();
   }
   else
   {
      dlib::matrix<double,1,1> intercept;
      intercept(0) = mean(y);

      if (continuous)
      {
         LinearLikelihood likelihood_fit(x, y);
         double SSE =  accu(square((y - x * intercept(0))));
         null_ll = -likelihood_fit(intercept) / (SSE/x.n_rows);
      }
      else
      {
         LogitLikelihood likelihood_fit(x, y);
         null_ll = -likelihood_fit(intercept);
      }
   }
   return null_ll;
}

// Likelihood-ratio test
double likelihoodRatioTest(Kmer& k, const double null_ll)
{
   double log_likelihood = k.log_likelihood();
   double lrt_p = 0;
   if (log_likelihood == 0 || null_ll == 0)
   {
      k.add_comment("zero-ll");
   }
   else
   {
      double lrt = 2*(log_likelihood - null_ll);
      if (lrt > 0)
      {
         lrt_p = normalPval(pow(lrt,0.5));
      }
   }
   return lrt_p;
}

// Returns p-value for a test statistic that is >0 and standard normally distributed
double normalPval(double testStatistic)
{
   double p_val = 0;
   if (testStatistic < 5)
   {
      boost::math::normal s;

      p_val = 2 * (1 - boost::math::cdf(s, testStatistic));
   }
   else
   {
      // For large z need to use a bound
      // See http://stats.stackexchange.com/questions/13690/how-to-compute-the-probability-associated-with-absurdly-large-z-scores
      //
      // Upper bound
      // S(z) <= phi(z)/z
      // cdf = 1-(0.5 * S(z))
      // At z = 5 correct to +/- 2.5%
#ifdef SEER_DEBUG
      std::cerr << "using erfc bound rather than 'exact' function\n";
#endif
      p_val = 2 * exp(-0.5*pow(testStatistic,2))*normalArea/testStatistic;
   }

   return p_val;
}
