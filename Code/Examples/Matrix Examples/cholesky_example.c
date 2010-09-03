#include <iostream>
#include <stdio.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;

int main(void){
  // Set the data for our matrix.
  
  double a[] = {2.0, 1.0, 0.0,
		1.0, 3.0, 2.0,
		0.0, 2.0, 4.0};


  // Now create a matrix structure using this data.
  gsl_matrix_view A = gsl_matrix_view_array(a, 3, 3);

  // Find the Cholesky Decomposition.
  gsl_linalg_cholesky_decomp(&A.matrix);

  /* NOTE:
     The Cholesky decomposition takes A and decomposes
     it according to L*L' = A where L is a lower triangular
     matrix.  Equivalently, U'*U = A where U is an upper
     triangular matrix and L = U'.  The gsl function above
     stores the upper _and_ lower triangular matrices U'
     and U in A.  Thus one needs to use the triangular
     multiplication rutines when multiplying by the output.
  */
  
  gsl_matrix_fprintf (stdout, &A.matrix, "%g"); cout<<"\n";

}
