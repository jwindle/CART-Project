#include <iostream>
#include <ctime>             // For time();
#include <gsl/gsl_rng.h>     // For gsl_rng
#include <gsl/gsl_randist.h> // For gsl_ran_gaussian()

using namespace std;

int main(void){
  // Generate a random number generator.
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
  // Seed the random number generator.
  // Time(NULL) returns seconds since Jan 1, 1970.
  gsl_rng_set(r, time(NULL));
  // Draw from a Gaussian(0,1)
  double draw = gsl_ran_gaussian(r, 1);
  cout<<"We drew "<<draw<<" from the Gaussian Urn.\n";
  gsl_rng_free(r);
  return 0;
}
