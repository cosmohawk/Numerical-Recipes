#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "myroots.h"

using namespace std;
//using namespace complex_literals;

/* rtbis begins with brackets for the root and then finds the root intself by
   bisection. Here we find all the roots of J_0 between 1 and 50. */

struct Funcd {
  // Structure containing the function and its derivative if known
  complex<double> operator() (const complex<double> z) {
    return pow(z,5.0) -1.0;
  }
  complex<double> df(const complex<double> z) {
    return 5.0*pow(z,4.0);
  }
};

int main(void) {

  complex<double> z1;
  z1.real(-4.0);
  z1.imag(4.0);
  cout << pow(z1,2.0) << endl;

  // starting point
  complex<double> z;
  complex<double> f;

  const int MAXIT = 100;
  const double tol = 1.0e-12;
  int i;

  Funcd fx;

  int m,n;
  for (m=0;m<2000;m++) {
    // Loop over real axis
    z1.real(m*0.0005);
    for (n=0;n<2000;n++) {
      // Loop over imaginary axis
      z1.imag(n*0.0005);
      
      // Check for NR convergence
      z = z1;
      for (i=0;i<MAXIT;i++){
	f = fx(z);
	if (abs(f) > pow(tol,2.0)) {
	  // Newton-Raphson step
	  z =  z - fx(z)/fx.df(z);
	} else {
	  // Print strting point if converged
	  cout << z1.real() << " " << z1.imag() << " " << 1 << " " <<i<< endl;
	  break;
	}
      }
      if (abs(f) > pow(tol,2.0))
	// Print staring point if not converged
	cout << z1.real() <<" "<< z1.imag()<<" "<<0<< " "<<100<<endl;
    }
  }
  
  return 0;
}
  
  
  

  
