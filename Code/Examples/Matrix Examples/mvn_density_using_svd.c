
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_blas.h>    // For matrix operations
#include <gsl/gsl_linalg.h>  // For SVD decomposition
#include <math.h>            // For sqrt()

using namespace std;

/*
  Note: Several other ways have come to mind about how one can
  evaluate the density of a normal random variable.  For instance, if
  one has the lower triangular cholesky decomposition of V, then he
  can do the following.  First, the determinant of V can be calculated
  as (the product of the diagonal entries of Chol)^2.  Second, if one
  has V^{-1}x = b, then x = Vb.  If we have the Cholesky decomposition
  we can solve x = Chol * Chol^T b by first solving x = Chol * c and
  then solving c = Chol^T * b, both by back substitution.  Then x^T
  V^{-1} x = x^T b.  Another idea: one can store both the variance
  matrix and the Cholesky decomposition in smaller data structures
  since we have symmetry and lower triangularity.
*/

/*
  Calculate the density of a multivariate normal random variable with
  mean zero and variance V at x.  The variance vector V can be
  decmoposed by SVD as UDU^T.  Hence the inverse of V can be
  decomposed as V^{-1} = UD^{-1/2}U^T.  If we let Sqrt = U D^{-1/2}
  then V^{-1} = Sqrt Sqrt^T.  Calculating x^T V^{-1} x is thus
  equivalent to (Sqrt^T x)^T (Sqrt^T x).  Furthermore, the determinant
  of V is identical to the determinant of D, which can easily be
  calculated by taking the product of the diagonal entries of D.
*/
int density(const gsl_matrix * Sqrt, const double det, const gsl_vector * x, 
	    double * result)
{
  // Return nonzero number for error.
  int error;

  // Get the length of the vector.
  unsigned int n = x->size;

  // If the number of columns in our matrix does not
  // match the number of elements in our vector then
  // we should return 0.
  if (Sqrt->size2 != n)
    return 1;

  /*
    Now use the SVD decomposition to calculate evaluate the density of
    a multivariate normal random variable at the point x.  As
    mentioned above The variance vector V can be decmoposed by SVD as
    UDU^T.  Hence the inverse of V can be decomposed as V^{-1} =
    UD^{-1/2}U^T.  If we let Sq = U D^{-1/2} then V^{-1} = Sq Sq^T.
    Calculating x^T V^{-1} x is thus equivalent to (Sq^T x)^T (Sq^T
    x). 
  */

  // Set up a temporary vector to work with.
  gsl_vector * temp = gsl_vector_alloc(n);

  // Calculate temp = 1.00 * Sqrt^T * x + 0.00 * temp.
  gsl_blas_dgemv (CblasTrans, 1.00, Sqrt, x, 0.00, temp);

  // Calcuate result = temp^T*temp = (Sqrt^T x)^T (Sqrt^T x).
  gsl_blas_ddot (temp, temp, result);

  // Calculate constant.
  double  K = exp( -0.5 * n * ( log(2*M_PI) + log(det) ) );

  // Calculate density at x.
  *result = K * exp(-0.5 * *result);
  
  return error;
}

int main(void){

  // Dimension of random variable.
  int n = 2;

  // Set the data for our covariance matrix.
  double cov_data[] = {1.0, 0.0,
		       0.0, 1.0};

  // Create a matrix structure using this data.
  gsl_matrix_view V = gsl_matrix_view_array(cov_data, n, n);

  /*
    We can calculate the inverse of V using the singular value
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
  double det = 1;

  /*
    Now construct the squareroot.  We do this by hand since I don't
    see an obvious BLAS routine to use.  Might as well get the
    determinant while we're at it.
  */
  // U*D^{1/2} = multiply each column u_j by S_{jj}^{1/2}.
  for(int j = 0; j < n; j++){
    // sqrtd = S_{jj}^{1/2} = ( S[j]^{1/2} )
    double sqrtd = sqrt( gsl_vector_get(S, j) );
    det = det * gsl_vector_get(S,j);
    // column v_j = u_j * S_{jj}^{1/2}.
    for(int i = 0; i< n; i++){
      // v_{ij} = u_{ij}*S_{jj}^{1/2}
      double v_ij = gsl_matrix_get(U, i, j) * sqrtd;
      gsl_matrix_set(SqrtV, i, j, v_ij);
    }
  }

  // Let's evaluate the density on some grid.

  /*
    Consider a one dimensional grid for the moment.  We have that x_0
    = x_start_value and x_K = x_end_value.  Thus there are K+1 data
    points on our grid with spacing given by delta = (x_end_value -
    x_start_value) / K.  Thus we need to create an array of K+1 data
    points and evaluate the grid at x_i = x_start_value + delta * i
    for i = 0, ..., K.
  */

  // The grid.
  int xgrid_size = 10;
  int ygrid_size = 10;
  // double * grid = new double[(xgrid_size+1) * (ygrid_size+1)];

  // The grid boundaries.
  double x_start_value = -1;
  double x_end_value = 1;
  double y_start_value = -1;
  double y_end_value = 1;

  // The grid size.
  double delta_x = (x_end_value - x_start_value) / xgrid_size;
  double delta_y = (y_end_value - y_start_value) / ygrid_size;

  // The point at which we evaluate the density.
  gsl_vector * input = gsl_vector_alloc(n);
  double output;

  for(int ix = 0; ix <= xgrid_size; ix++){
    for(int iy = 0; iy <= ygrid_size; iy++){
      double x = x_start_value + delta_x * ix;
      double y = y_start_value + delta_y * iy;
      gsl_vector_set(input, 0, x);
      gsl_vector_set(input, 1, y);
      density(SqrtV, det, input, &output);
      cout<<x<<" "<<y<<" "<<output<<"\n";
    }
  }
 
}
