/*
 * File: panglossStruct.cpp
 *
 * Implements metric MDS (multi-dimensional scaling)
 * for pangloss kmers
 */

#include "pangloss.hpp"

arma::mat metricMDS(const arma::mat& populationMatrix, const int dimensions)
{
   /*
    * Metric MDS
    *
    * 1) P^2 -> matrix with elements which are distances squared
    * 2) J = I - n^-1(II') - II' is a square matrix of ones
    * 3) B = -0.5JP^2J
    * 4) Decompose B into eigenvalues
    * 5) MDS components = eigenvectors * eigenvalues
    */
   const unsigned int matSize = populationMatrix.n_rows;

   // Step 1)
   arma::mat P = arma::square(dissimiliarityMatrix(populationMatrix));

   // Step 2)
   arma::mat J = arma::eye<arma::mat>(matSize, matSize)
      - 1/matSize*arma::ones<arma::mat>(matSize, matSize);

   // Step 3)
   arma::mat B = -0.5 * J * P * J;

   // Step 4)
   arma::vec eigval;
   arma::mat eigvec;

   arma::eig_sym(eigval, eigvec, B);

   // Step 5)
   arma::mat mds = eigvec * diagmat(eigval);

   // Sort eigenvalues
   arma::uvec sort_idx;
   sort_idx = sort_index(square(eigval), "descend");

   arma::mat mdsCols(matSize, dimensions);
   for (unsigned int i = 0; i < dimensions; ++i)
   {
      mdsCols.col(i) = mds.col(sort_idx(i));
   }

   return mdsCols;
}

// Distance between all rows. 0/1 elements only
arma::mat dissimiliarityMatrix(const arma::mat& inMat)
{
   const unsigned int matSize = inMat.n_rows;
   arma::mat dist(matSize, matSize);

   // Loop through upper triangle
   for (unsigned int i = 0; i < matSize; ++i)
   {
      arma::vec ref_row = inMat.row(i).t();
      for (unsigned int j = i; j < matSize; j++)
      {
         if (i == j)
         {
            dist(i, j) = 0;
         }
         else
         {
            dist(i, j) = accu(square(ref_row * inMat.row(j))); // Distance function
            dist(j, i) = dist(i, j); // Set symmetric elements
         }
      }
   }

   return dist;
}
