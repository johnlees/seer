/*
 * File: kmdsStruct.cpp
 *
 * Implements metric MDS (multi-dimensional scaling)
 * for kmds kmers
 */

#include "kmds.hpp"

arma::mat metricMDS(const arma::mat& populationMatrix, const int dimensions, const unsigned int threads, const std::string& distances_file)
{
   /*
    * Metric MDS
    *
    * 1) P^2 -> matrix with elements which are distances squared
    * 2) J = I - n^-1(II') - II' is a square matrix of ones
    * 3) B = -0.5JP^2J
    * 4) Decompose B into eigenvalues
    * 5) MDS components = eigenvectors * eigenvalues
    * 6) Normalise components
    */
   const unsigned int matSize = populationMatrix.n_rows;

   // Step 1)
   arma::mat P = arma::square(dissimiliarityMatrix(populationMatrix, threads));

   // If supplied as an optional parameter, write distance matrix to file
   if (!distances_file.empty())
   {
      writeDistances(distances_file, P);
   }

   // Step 2)
   arma::mat J = arma::eye<arma::mat>(matSize, matSize)
      - 1/double(matSize)*arma::ones<arma::mat>(matSize, matSize);

   // Step 3)
   arma::mat B = -0.5 * J * P * J;

   // Step 4)
   arma::vec eigval;
   arma::mat eigvec;

   arma::eig_sym(eigval, eigvec, B);

   // Step 5)
   // Eigenvalues returned are sorted ascending, so want to reverse order
   arma::mat mds = fliplr(eigvec * diagmat(sqrt(eigval)));

   // Step 6)
   // All values will lie in the interval [-1,1]
   arma::mat norm_mds(matSize, dimensions);
   for (int i = 0; i < dimensions; ++i)
   {
      if (pow(max(mds.col(i)), 2) > pow(min(mds.col(i)), 2))
      {
         norm_mds.col(i) = mds.col(i)/max(mds.col(i));
      }
      else
      {
         norm_mds.col(i) = mds.col(i)/fabs(min(mds.col(i)));
      }
   }

   return norm_mds;
}

// Distance between all rows. 0/1 elements only
arma::mat dissimiliarityMatrix(const arma::mat& inMat, const unsigned int threads)
{
   const unsigned int matSize = inMat.n_rows;
   arma::mat dist = arma::zeros<arma::mat>(matSize, matSize);

   // Time parallelisation
   std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

   // Create queue for distance calculations
   std::queue<std::future<std::vector<DistanceElement>>> distance_calculations;

   // Set up number of jobs per thread
   const unsigned long int num_calculations = 0.5*matSize*(matSize - 1);
   const unsigned long int calc_per_thread = (unsigned long int)num_calculations / threads;
   const unsigned int num_big_threads = num_calculations % threads;

   // Keep track of row and column outside of outer loop
   unsigned int row = 0;
   unsigned int col = 0;
   for (unsigned int thread_idx = 0; thread_idx < threads; ++thread_idx) // Loop over threads
   {
      // First 'big' threads have an extra job
      unsigned long int thread_jobs;
      if (distance_calculations.size() < num_big_threads)
      {
         thread_jobs = calc_per_thread + 1;
      }
      else
      {
         thread_jobs = calc_per_thread;
      }

      // Each thread is assigned elements to calculated scanning left to right,
      // top to bottom on the upper triangle of the matrix
      std::vector<DistanceElement> thread_elements;
      thread_elements.reserve(thread_jobs);
      for (unsigned int element = 0; element < calc_per_thread; ++element)
      {
         DistanceElement d;
         d.row = row;
         d.col = ++col;

         thread_elements.push_back(d);

         if (col + 1 == matSize)
         {
            ++row;
            col = row;
         }
      }

      // Set the thread off
      distance_calculations.push(std::async(std::launch::async, threadDistance, thread_elements, std::cref(inMat)));
   }

   // Wait for thread results to return, and put the results in the distance
   // matrix
   while(!distance_calculations.empty())
   {
      std::vector<DistanceElement> distance_elements = distance_calculations.front().get();
      distance_calculations.pop();

      for(std::vector<DistanceElement>::iterator it = distance_elements.begin() ; it < distance_elements.end(); ++it)
      {
         dist(it->row, it->col) = it->distance;
         dist(it->col, it->row) = it->distance;
      }
   }

   // Print time taken
   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
   std::chrono::duration<double> diff = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
   std::cerr << "Distance matrix calculated in: " << diff.count() << " s\n";

   return dist;
}

std::vector<DistanceElement> threadDistance(std::vector<DistanceElement> element_list, const arma::mat& rectangular_matrix)
{
   for(std::vector<DistanceElement>::iterator it = element_list.begin() ; it < element_list.end(); ++it)
   {
      it->distance = distanceFunction(rectangular_matrix.row(it->row), rectangular_matrix.row(it->col));
   }

   return element_list;
}

double distanceFunction(const arma::rowvec& vec_1, const arma::rowvec& vec_2)
{
   // Slow
   // return accu(abs(vec_1 - vec_2));
   return dot(vec_1 - vec_2, vec_1 - vec_2);
}
