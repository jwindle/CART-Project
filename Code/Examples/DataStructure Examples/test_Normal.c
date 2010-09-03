#include <iostream>
#include <stdio.h>
#include <gsl/gsl_blas.h>    // For matrix operations
#include <gsl/gsl_linalg.h>  // For Cholesky Decomp
#include <ctime>             // For time()
#include <gsl/gsl_rng.h>     // For gsl_rng
#include <gsl/gsl_randist.h> // For gsl_ran_gaussian()
#include <gsl/gsl_sf_exp.h>  // For gsl_sf_exp()
#include "Normal.h"          // For Normal

using namespace std;

double SqrExpCovFnc(double ti, double tj, double k1, double k2, double k3)
{
  return k1*gsl_sf_exp(-0.5 * (ti-tj)*(ti-tj) / k2) + k3*(ti==tj);
}

int main(void){

  /*
    Generate a random number generator.  There are several different
    types of random number generators one may use.  I simply chose
    gsl_rng_taus since that is the example gsl uses.
  */
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

  // We want to produce a new sequence of random
  // numbers each time we run this program.  Thus
  // we need to seed the random number generator.

  // Seed the random number generator.
  // Time(NULL) returns seconds since Jan 1, 1970.
  gsl_rng_set(r, time(NULL));

  // The number of time steps we will produce.
  int N = 3;

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

  // Create a convariance matrix using this data.
  gsl_matrix_view V = gsl_matrix_view_array(cov_data, N, N);

  // Set the mean to zero.
  gsl_vector * m;
  m = gsl_vector_alloc(N);
  gsl_vector_set_zero(m);

  cout<<"I got here.\n";

  // Now create our normal random variable.
  Normal * normal;
  normal = new Normal(m, &V.matrix);

  // Now let's make several draws from a normal
  // random variable with covariance matrix V.

  gsl_vector * draw;
  draw = gsl_vector_alloc(N);

  for(int j = 0; j< 100000; j++){
    normal->sample(r, draw);
    // See what we get.
    for(int i = 0; i<N; i++)
  	cout<<gsl_vector_get(draw, i)<<" ";
    cout<<"\n";
  }

  return 0;
}
