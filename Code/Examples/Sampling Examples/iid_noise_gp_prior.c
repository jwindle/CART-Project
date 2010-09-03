#include <iostream>
#include <stdio.h>
#include <gsl/gsl_blas.h>    // For matrix operations
#include <gsl/gsl_linalg.h>  // For Cholesky Decomp
#include <ctime>             // For time()
#include <gsl/gsl_rng.h>     // For gsl_rng
#include <gsl/gsl_randist.h> // For gsl_ran_gaussian()
#include <gsl/gsl_sf_exp.h>  // For gsl_sf_exp()

using namespace std;

double SqrExpCovFnc(double ti, double tj, double k1, double k2, double k3)
{
  return k1*gsl_sf_exp(-0.5 * (ti-tj)*(ti-tj) / k2) + k3*(ti==tj);
}

int DrawNormal(const gsl_matrix * Chol, const gsl_rng * r, gsl_vector * draw)
{
  // Record the status of our procedure.
  int error = 0;

  // Get the length of the vector.
  int n = draw->size;

  // If the number of columns in our matrix does not
  // match the number of elements in our vector then
  // we should return 0.
  if (Chol->size2 != draw->size)
    return 1;

  // We take n independent draws from a Gaussian.
  double random_num;
  for(int i = 0; i < n; i++){
    // Draw from a Gaussian(0,1)
    random_num = gsl_ran_gaussian(r, 1);
    gsl_vector_set(draw, i, random_num);
  }

  // Now use the Cholesky decomposition to draw from a 
  // multivariate normal with covariance matrix Chol*Chol^T.

  // Calculate draw = Chol * iid.
  // We need to be careful here since Chol is supposed to be
  // lower triangular.

  // For instance, this is wrong, it is full matrix multiplication.
  // It would be correct if gsl did what we thought it should.
  // error = gsl_blas_dgemv(CblasTrans, 1.0, Chol, iid, 0.0, draw);

  // This is right, it multiplies our vector only using the lower
  // triangular entries of Chol.  It is bascially x = Chol * x.
  gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit, Chol, draw);

  return error;
}

int main(void){

  // Generate a random number generator.
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

  // We want to produce a new sequence of random
  // numbers each time we run this program.  Thus
  // we need to seed the random number generator.

  // Seed the random number generator.
  // Time(NULL) returns seconds since Jan 1, 1970.
  gsl_rng_set(r, time(NULL));

  // The number of time steps we will produce.
  int N = 10;

  // The parameters for our covariance matrix.
  double k1 = 0.5;
  double k2 = 1.0;
  double k3 = 0.5;

  // Set aside space for our covariance matrix.
  double * cov_data = new double[N*N];

  // Set the data for our covariance matrix.
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      cov_data[i*N+j] = SqrExpCovFnc(i, j, k1, k2, k3);

  // Create a matrix structure using this data.
  gsl_matrix_view V = gsl_matrix_view_array(cov_data, N, N);

  // To check our matrix:
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++)
      cout<<gsl_matrix_get(&V.matrix, i, j)<<" ";
    cout<<"\n";
  }
  /*
    Suppose X is a multivariate normal RV with covariance matrix I.
    Then we can construct a random variable Y with covariance matrix V
    by using X and the Cholesky decomposition of V.  In particular, if
    AA^T = V, then E[(AX)(AX)^T] = A E[XX^T] A^T = V.

    Note that the gsl implementation of the Cholesky decomposition
    differs from the Octave implementation.  In octave A = chol(V)
    solves A^T A = V where A is uppper triangular .  With gsl the
    Cholesky decomposition A contains both the upper and lower
    triangular matrices A and A^T.  Thus we must use a triangular
    matrix multiplcation routine when generating our multivariate
    normal draw.
  */

  // We need the Cholesky decomposition of V.
  gsl_linalg_cholesky_decomp(&V.matrix);

  /*
    The decomposition is stored in V.matrix where the lower triangular
    entries of A are stored in the lower triangular entries of V and
    the upper triangular entries of A^T are stored in the upper
    triangular entries of V.
  */

  // Now let's make several draws from a normal
  // random variable with covariance matrix V.

  gsl_vector * draw;
  draw = gsl_vector_alloc(N);
  
  for(int j = 0; j< 100000; j++){
    DrawNormal(&V.matrix, r, draw);
    // See what we get.
    for(int i = 0; i<N; i++)
  	cout<<gsl_vector_get(draw, i)<<" ";
    cout<<"\n";
  }

  /*
    You can check this data by piping the output of this program into
    a text file, for instance normal.data, and then calculating the
    sample variance of the data.  It should look like the covariance
    matrix V.
  */

}
