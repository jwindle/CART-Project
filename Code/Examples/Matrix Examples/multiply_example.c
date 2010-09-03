#include <iostream>
#include <stdio.h>
#include <gsl/gsl_blas.h>

using namespace std;

int main(void){

  // Set the data for our matrix.
  double A_data[] = {1.0, 1.0, 0.0,
		1.0, 2.0, 1.0,
		0.0, 1.0, 4.0};

  // Set the data for our vector.
  double x_data[] = {1.0, 3.0, 2.0};

  // Now create a matrix structure using this data.
  gsl_matrix_view A = gsl_matrix_view_array(A_data, 3, 3);

  // Create a vector structure using this data.
  gsl_vector_view x = gsl_vector_view_array(x_data, 3);

  // Allocate space for y in y = Ax.
  gsl_vector * y = gsl_vector_alloc(3);

  // Calculate y = Ax.
  gsl_blas_dgemv (CblasNoTrans, 1.0, &A.matrix, &x.vector, 0.0, y);

  // Print the output.
  gsl_matrix_fprintf (stdout, &A.matrix, "%g"); cout<<"\n";
  gsl_vector_fprintf (stdout, &x.vector, "%g"); cout<<"\n";
  gsl_vector_fprintf (stdout, y, "%g");

  // Free up the memory we dynamically allocated.
  gsl_vector_free(y);
}
