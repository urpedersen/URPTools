#include "randomGaussian.h"

#include <cstdlib>
#include <cmath>

/**
 * Generate a random number using the polar Box–Muller transform method (aka Marsaglia polar method).
 */
double randomGaussian() {
  double u,v,s = 2.0;
  do {
    u = 2.0*(double)rand()/(double)RAND_MAX-1.0;
    v = 2.0*(double)rand()/(double)RAND_MAX-1.0;
	s = u*u + v*v;
  } while( s >= 1.0 || s == 0.0 );
  return u*sqrt(-2.0*log(s)/s);  
  // Note another random number is avalible in v*sqrt(-2.0*log(s)/s)
  
  // Below is an implimentation using a sum of 12 random numbers
  //double out = 0.0;
  //for(unsigned i=0;i<12;i++) out+=((double)rand()/(double)RAND_MAX-0.5); // Normal distribution with sigma=1
  //return out;
}
