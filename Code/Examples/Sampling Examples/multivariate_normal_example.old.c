#include <iostream>
#include <stdio.h>
#include <gsl/gsl_blas.h>    // For matrix operations
#include <gsl/gsl_linalg.h>  // For Cholesky Decomp
#include <ctime>             // For time()
#include <gsl/gsl_rng.h>     // For gsl_rng
#include <gsl/gsl_randist.h> // For gsl_ran_gaussian()

using namespace std;

int main(void){

  // Generate a random number generator.
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

  // We want to produce a new sequence of random
  // numbers each time we run this program.  Thus
  // we need to seed the random number generator.

  // Seed the random number generator.
  // Time(NULL) returns seconds since Jan 1, 1970.
  gsl_rng_set(r, time(NULL));

  // Set the data for our covariance matrix.
  double cov_data[] = {1.0, 1.0, 0.0,
		       1.0, 2.0, 1.0,
		       0.0, 1.0, 4.0};

  // Create a matrix structure using this data.
  gsl_matrix_view C = gsl_matrix_view_array(cov_data, 3, 3);

  /*
    Suppose X is a multivariate normal RV with covariance
    matrix I.  Then we can construct a random variable Y
    with covariance matrix C by using X and the Cholesky
    decomposition of C.  In particular, if AA^T = C, then
    E[(AX)(AX)^T] = A E[XX^T] A^T = C.
  */

  // We need the Cholesky decomposition of C.
  gsl_linalg_cholesky_decomp(&C.matrix);

  // We need to draw a 3 dimensional RV.
  double draw_data[3];

  // We take 3 independent draws from a Gaussian.
  for(int i = 0; i < 3; i++){
    // Draw from a Gaussian(0,1)
    draw_data[i] = gsl_ran_gaussian(r, 1);
  }

  // Load this data into a vector structure.
  gsl_vector_view iid_draw = gsl_vector_view_array(draw_data, 3);

  // Now use the Cholesky decomposition to draw
  // from a multivariate normal with covariance matrix C.

  // Allocate space for draw in draw = A iid_draw.
  gsl_vector * draw = gsl_vector_alloc(3);

  // Calculate draw = A iid_draw.
  gsl_blas_dgemv(CblasNoTrans, 1.0, &C.matrix, &iid_draw.vector, 0.0, draw);

  // See what we get.
  for(int i = 0; i<3; i++)
    cout<<gsl_vector_get(draw, i)<<" ";

  cout<<"\n";

}
