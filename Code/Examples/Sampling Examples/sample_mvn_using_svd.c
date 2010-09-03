
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_blas.h>    // For matrix operations
#include <gsl/gsl_linalg.h>  // For SVD decomposition
#include <ctime>             // For time()
#include <gsl/gsl_rng.h>     // For gsl_rng
#include <gsl/gsl_randist.h> // For gsl_ran_gaussian()
#include <math.h>            // For sqrt()

using namespace std;

/*
  Draw a normal random variable with variance V = Sqrt*Sqrt^T, where
  Sqrt is the ``squareroot'' of V as constructed using the SVD.  In
  particular, when V = U*D*U^T, then Sqrt = U*D^{1/2}.
*/
int DrawNormal(const gsl_matrix * Sqrt, const gsl_rng * r, gsl_vector * draw)
{
  // Record the status of our procedure.
  int error = 0;

  // Get the length of the vector.
  int n = draw->size;

  // If the number of columns in our matrix does not
  // match the number of elements in our vector then
  // we should return 0.
  if (Sqrt->size2 != draw->size)
    return 1;

  // Get space to draw a collection of iid random variable.
  gsl_vector * iid;
  iid = gsl_vector_alloc(n);

  // We take n independent draws from a Gaussian.
  double random_num;
  for(int i = 0; i < n; i++){
    // Draw from a Gaussian(0,1)
    random_num = gsl_ran_gaussian(r, 1);
    gsl_vector_set(iid, i, random_num);
  }

  /*
    Now use the SVD decomposition to draw from a multivariate normal
    random variable.  We can decompose the variance V = U*D*U^T using
    the SVD where U is a unitary matrix and D is diagonal.  If we let
    Sqrt = U*D^{1/2} then V = Sqrt*Sqrt^T.  We can use this matrix to
    construct a random variable with variance V.  In particular, if X
    is an n-dimensional standard normlal random variable then Y =
    Sqrt*X is a random variable with variance E[YY^T] = Sqrt*Sqrt^T =
    V.
  */

  // Calcuate draw = 1.00 * Sqrt * iid + 0.00 * draw.
  gsl_blas_dgemv (CblasNoTrans, 1.00, Sqrt, iid, 0.00, draw);

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
  
  // Dimension of random variable.
  int n = 3;

  // Set the data for our covariance matrix.
  double cov_data[] = {2.0, 1.0, 0.0,
		       1.0, 3.0, 2.0,
		       0.0, 2.0, 4.0};

  // Create a matrix structure using this data.
  gsl_matrix_view V = gsl_matrix_view_array(cov_data, n, n);

  /*
    Suppose X is a multivariate normal RV with covariance matrix I.
    Then we can construct a random variable Y with covariance matrix V
    by using X and the Singular Value Decomposition of V.  In
    particular, if V = UDU^T, then we can define a ``sqareroot'' of V
    by Sq = U D^{1/2} so that V = Sq*Sq^T.  The random variable Y =
    Sq*X then has variance E[Y*Y*T] = Sq*Sq^T = V.

    We can also calculate the inverse of V using the singular value
    decomposition.  Again, suppose V = UDU^T.  Then V^{-1} =
    UD^{-1}U^T.  If we let now let Sq = U D^{-1/2}, then we have that
    V^{-1} = Sq*Sq^T.  This is useful when we want to evaluate the
    density of a multivariate normal random variable.
  */

  /*
    We need the singular value decomposition of V.  In this case, the
    gsl_linalg_SV_decomp function factorizes V into the singular value
    decomposition A = U S U^T. On output the matrix V is replaced by
    U. The diagonal elements of the singular value matrix S are stored
    in the vector S a workspace of length N = (# cols of V) is
    required in work.
  */
  gsl_matrix * U = gsl_matrix_alloc(n, n);
  gsl_vector * S = gsl_vector_alloc(n);
  gsl_vector * work = gsl_vector_alloc(n);

  gsl_linalg_SV_decomp(&V.matrix, U, S, work);

  /* 
     Now lets consruct our ``squareroot.'' I put squareroot in
     quotations because there aer several ways one can construct a
     squareroot.
  */
  gsl_matrix * SqrtV = gsl_matrix_alloc(n, n);

  /*
    Now construct the squareroot.  We do this by hand since I don't
    see an obvious BLAS routine to use.
  */
  // U*D^{1/2} = multiply each column u_j by S_{jj}^{1/2}.
  for(int j = 0; j < n; j++){
    // sqrtd = S_{jj}^{1/2} = ( S[j]^{1/2} )
    double sqrtd = sqrt( gsl_vector_get(S, j) );
    // column v_j = u_j * S_{jj}^{1/2}.
    for(int i = 0; i< n; i++){
      // v_{ij} = u_{ij}*S_{jj}^{1/2}
      double v_ij = gsl_matrix_get(U, i, j) * sqrtd;
      gsl_matrix_set(SqrtV, i, j, v_ij);
    }
  }

  // Now let's make several draws from a normal
  // random variable with covariance matrix V.

  gsl_vector * draw;
  draw = gsl_vector_alloc(n);
  
  for(int j = 0; j<10000; j++){
    DrawNormal(SqrtV, r, draw);
    // See what we get.
    for(int i = 0; i<n; i++)
  	cout<<gsl_vector_get(draw, i)<<" ";
    cout<<"\n";
  }

  /*
    You can check this data by piping the output of this program into
    a text file, for instance normal.data, and then calculating the
    sample variance of the data.  It should look like the covariance
    matrix V.

    The R code is
    y = read.table('normal.data');
    cov(y);
  */

}
