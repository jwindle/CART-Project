#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_blas.h>    // For matrix operations
#include <gsl/gsl_linalg.h>  // For Cholesky Decomp
#include <ctime>             // For time()
#include <gsl/gsl_rng.h>     // For gsl_rng
#include <gsl/gsl_randist.h> // For gsl_ran_gaussian()
#include <gsl/gsl_sf_exp.h>  // For gsl_sf_exp()

using namespace std;

/*
  gsl_matrix_printf prints a matrix as a column vector.  This function
  prints a matrix in block form.
*/
void pretty_print(const gsl_matrix * M)
{
  // Get the dimension of the matrix.
  int rows = M->size1;
  int cols = M->size2;
  // Now print out the data in a square format.
  for(int i = 0; i < rows; i++){
    for(int j = 0; j < cols; j++){
      cout<<gsl_matrix_get(M, i, j)<<" ";
    }
    cout<<"\n";
  } 
}

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

  /*
    Suppose X is a multivariate normal RV with covariance matrix I.
    Then we can construct a random variable Y with covariance matrix K
    by using X and the Cholesky decomposition of K.  In particular, if
    AA^T = K, then E[(AX)(AX)^T] = A E[XX^T] A^T = K.

    Note that the gsl implementation of the Cholesky decomposition
    differs from the Octave implementation.  In octave A = chol(K)
    solves A^T A = K where A is uppper triangular .  With gsl the
    Cholesky decomposition A contains both the upper and lower
    triangular matrices A and A^T.  Thus we must use a triangular
    matrix multiplcation routine when generating our multivariate
    normal draw.
  */

  // This is right, it multiplies our vector only using the lower
  // triangular entries of Chol.  It is bascially x = Chol * x.
  gsl_blas_dtrmv (CblasLower, CblasNoTrans, CblasNonUnit, Chol, draw);

  return error;
}

// x' A x
double SymQuadraticForm(const gsl_matrix * A, const gsl_vector * x)
{
  double result;
  gsl_vector * temp_vec = gsl_vector_alloc(A->size2);
  gsl_blas_dsymv(CblasUpper, 1.0, A, x, 0.0, temp_vec);
  gsl_blas_ddot(x, temp_vec, &result);
  delete temp_vec;
  return result;
}

void GibbsSampling(const gsl_vector * y,
		  const gsl_matrix * KInvQ,
		  const gsl_matrix * CholC,
		  const gsl_matrix * InvA,
		  const gsl_rng * r,
		  const double n_prior,
		  const double d_prior,
		  gsl_vector * phi_data)
{
  // The size of our sample.
  int num_samples = phi_data->size;

  // The length of our time series.
  int T = y->size;
  
  // A few preliminary calculations.
  
  // The mean of the conditional density p(f | y,sig).
  gsl_vector * m = gsl_vector_alloc(T);
  gsl_blas_dgemv(CblasNoTrans, 1.0, KInvQ, y, 0.0, m);

  /*
    We must start our Gibbs sampling with some value for phi.  We will
    simply pull from the prior distribution to start.
  */
  double draw_gamma = gsl_ran_gamma(r, n_prior/2, d_prior/2);
  gsl_vector_set(phi_data, 0, draw_gamma);

  // Set aside space to draw a normal random variable.
  gsl_vector * draw = gsl_vector_alloc(T);

  // Set aside space for x = (y,f);
  gsl_vector * x = gsl_vector_alloc(2*T);
  for(int i = 0; i < T; i++){
    gsl_vector_set(x, i, gsl_vector_get(y, i) );
  }

  // We know beforehand n_post_joint.
  double n_post_joint = n_prior + 2 * T;

  /* THE ACTUAL SAMPLING TAKES PLACE HERE */
  for(int i = 1; i < num_samples; i++){
    /* DRAW NORMAL */
    // We want to draw f|phi ~ N(m, C/phi).
    // The derivation can be found in J. Scott's notes.
    
    // We break this into several steps.
    // First, f ~ N(0,C).
    DrawNormal(CholC, r, draw);
    // Second, f = f/sqrt(phi) ~ N(0, C/phi).
    double s = 1.0/sqrt(gsl_vector_get(phi_data,i-1));
    // cout<<1.0/s<<"\t";
    gsl_vector_scale(draw, s);
    // Third, f = f + m ~ N(m, C/phi).
    gsl_vector_add(draw, m);
    
    /* DRAW GAMMA */
    // Now we draw from phi | (y,f).
    
    // Let x = (y,f).
    for(int k = 0; k < T; k++){
      gsl_vector_set(x, k+T, gsl_vector_get(draw, k) );
    }
    /*
    for(int k = 0; k< 2*T; k++)
      cout<<gsl_vector_get(x,k)<<"\t";
    cout<<"\n";
    */
    /*
    for(int k=0; k < 2*T; k++){
      cout<<gsl_vector_get(x,k)<<"\t";
    }
    cout<<"\n";
    */
    double d_post_joint = SymQuadraticForm(InvA, x);

    double result;
    gsl_blas_ddot(y, draw, &result);
    // cout<<result<<"\t";
    // cout<<d_post_joint<<"\t";
       
    d_post_joint = d_post_joint + d_prior;
    //  cout<<d_post_joint<<"\t";

  
    // Now draw phi | (y,f).
    draw_gamma = gsl_ran_gamma(r, n_post_joint/2.0, 2.0/d_post_joint);
    // cout<<sqrt(draw_gamma)<<"\t";
    gsl_vector_set(phi_data, i, draw_gamma);

    // cout<<"\n"; 
  }

}

void ReadDataFile(const char * filename, double * & array, int & array_size)
{
    ifstream instream;
    
    instream.open(filename);
    if (!instream) {
        cout << "Unable to open file.\n";
        exit(1); // Stop things immediately.
    }

    double buffer[10000];
    array_size = 0;

    double value;

    while ( (instream >> value) && (array_size < 10000) ) {
      buffer[array_size] = value;
      array_size += 1;
    }

    // if (array != NULL) delete array;

    array = new double[array_size];

    for(int i = 0; i < array_size; i++)
      array[i] = buffer[i];
    
    instream.close();
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
  // int T = 2;

  // The parameters for our covariance matrix.
  /*
  double k1 = 0.5;
  double k2 = 1.0;
  double k3 = 0.5;
  */

  // The prior parameters.
  double n_prior = 3.0;
  double d_prior = 2.0;

  double * y_data;
  int T;

  ReadDataFile("y.data", y_data, T);
  //cout<<T<<"\t"<<y_data[0]<<"\t"<<y_data[1]<<"\n";

  /*
  gsl_vector * y = gsl_vector_alloc(T);
  gsl_vector_set(y, 0, -1.569684);
  gsl_vector_set(y, 1, -1.705319);
  */

  gsl_vector_view y_view = gsl_vector_view_array(y_data, T);

  /* CREATE VARIOUS MATRICES */

  /* THIS TAKES A FEW LINES OF CODE... */

  /*
    To remind the reader.  We are working with the following
    probability model: 
      y|f,phi ~ N(f, I/phi); 
      f|phi   ~ N(0, K/phi);
      phi     ~ Ga(n/2,d/2); 

    We assume that we have observed one data point y.  We want to
    calculate the posterior distribution for phi using Gibbs sampling.
    More precisely, we will calculate the joint distribution of
    (f,phi) | y using Gibbs samplingg and then marginalize to find the
    posterior distribution for phi | y.

    To find the joint distribution we draw repeatedly from the
    conditional ditributions p(f | y,\phi) and p(\phi | y,f).  In the
    limit this is as if we draw from p(phi, f | y).
  */

  // We start with the covariance matrix K for which
  // f | sig ~ N(0, sig^2 K).

  /*
  // Set aside space for our covariance matrix.
  double * cov_data = new double[T*T];

  // Set the data for our covariance matrix.
  for(int i = 0; i < T; i++)
    for(int j = 0; j < T; j++)
      cov_data[i*T+j] = SqrExpCovFnc(i, j, k1, k2, k3);
  */
  double * cov_data;
  int cov_data_size;

  ReadDataFile("K.data", cov_data, cov_data_size);

  if ( (T * T) != cov_data_size ){
    cout<<"Observation vector does not match size of covariance matrix.\n";
    exit(1);
  }

  // Create a matrix K using this data.
  gsl_matrix_view K = gsl_matrix_view_array(cov_data, T, T);

  // pretty_print(&K.matrix);

  // To check our matrix:
  /*
  for(int i = 0; i < T; i++){
    for(int j = 0; j < T; j++)
      cout<<gsl_matrix_get(&K.matrix, i, j)<<" ";
    cout<<"\n";
  }
  */

  // The identity matrix.
  gsl_matrix * Id = gsl_matrix_alloc(T, T);
  gsl_matrix_set_identity(Id);

  // The variance of y ~ N(0, Q = I + K).
  gsl_matrix * Q = gsl_matrix_alloc(T, T);
  // Initially, set Q = I + K;
  gsl_matrix_memcpy(Q, &K.matrix); // Q = K;
  gsl_matrix_add (Q, Id);          // Q = Q + I;

  // We need to calculate KInvQ.
  gsl_linalg_cholesky_decomp(Q);
  gsl_matrix * KInvQ = gsl_matrix_alloc(T,T);
  // KInvQ = K;
  gsl_matrix_memcpy(KInvQ, &K.matrix);
  // KInvQ = KInvQ * Chol(Q)^{-T}.
  gsl_blas_dtrsm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 
		 1.0, Q, KInvQ);
  // KInvQ = KInvQ * Chol(Q)^{-1}.
  gsl_blas_dtrsm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 
		 1.0, Q, KInvQ);

  // We need to calculate C = K - KInvQ K.
  gsl_matrix * C = gsl_matrix_alloc(T, T);
  gsl_matrix_memcpy(C, &K.matrix);
  // C = -1.0 * KInvQ * K + 1.0 * C.
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, KInvQ, &K.matrix, 1.0, C);

  // We need the Cholesky decomposition of C to sample a normal.

  /*
    The decomposition is stored in C where the lower triangular
    entries of C (original) are stored in the lower triangular entries
    of C (updated) and the upper triangular entries of Chol(C)^T are
    stored in the upper triangular entries of C (updated).

    I am familiar with the Cholesky decomposition D = Chol(C) where D
    is a lower triangular matrix and DD^T = C.  We want the lower
    triangular portion of D so that D*X, where X is a standard normal,
    has variance DD^T = C.  We have carefully chosen the correct
    version of BLAS matrix multiplication for this to work in
    DrawNormal().
  */
  gsl_linalg_cholesky_decomp(C);
  
  // Now we can find the inverse of C.
  // This will be useful to calculate InvA.
  gsl_matrix * InvC = gsl_matrix_alloc(T, T);
  // InvC = Id.
  gsl_matrix_memcpy(InvC, Id);
  // InvC = InvC * Chol(C)^{-T}
  gsl_blas_dtrsm(CblasRight, CblasLower, CblasTrans, CblasNonUnit, 
		 1.0, C, InvC);
  // InvC = InvC * Chol(C)^{-1}.
  gsl_blas_dtrsm(CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, 
		 1.0, C, InvC);

  /*
    The joint distribution of (y,f) | sig is normal with mean zero and
    variance A given by A = [Q, K; K, K]; In our calculations we do
    not need A, but rather A^{-1}.  If you refer to page 33 of my
    (June 30) notes you will see that InvA = [I, -I; -I, C^{-1}] where
    C = K - K InvQ K.  This is the same quantity we need to calculate
    the variance of f | y, sig.
  */ 

  // Notice we can simplify the calculation x' InvA x.  See notes.
   
  gsl_matrix * InvA = gsl_matrix_alloc(2*T, 2*T);

  for(int i = 0; i < T; i++){
    for(int j = 0; j < T; j++){
      double Id_ij = gsl_matrix_get(Id, i, j);
      double InvC_ij = gsl_matrix_get(InvC, i, j);
      gsl_matrix_set (InvA, i, j, Id_ij);
      gsl_matrix_set (InvA, i+T, j, -1.0*Id_ij);
      gsl_matrix_set (InvA, i, j+T, -1.0*Id_ij);
      gsl_matrix_set (InvA, i+T, j+T, InvC_ij);
    }
  }

  // Check that we did things right.
  //pretty_print(&K.matrix);
  //pretty_print(Q);
  //pretty_print(InvA);
  //pretty_print(InvC);  

  /* DONE CREATING NECESSARY MATRICES */

  int num_samples = 5000;
  gsl_vector * phi_data = gsl_vector_alloc(num_samples);
  
  GibbsSampling(&y_view.vector,
		KInvQ,
		C,
		InvA,
	        r,
		n_prior,
		d_prior,
	        phi_data);

  for(int i=0; i < num_samples; i++){
    cout<<gsl_vector_get(phi_data,i)<<"\n";
  }

}
