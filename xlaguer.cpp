#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include "../code/nr3.h"
//#include "../code/roots_poly.h"
#include "myroots_poly.h"


using namespace std;

/* laguer finds the roots of a polynomial with complex coefficients. The 
   polynomial in this case is: 

   F(x) = x^4 - (1+2i)x^2 +2i.

   The four roots of the polynomial are x=1.0, -1.0, 1+i, and -(1+i). */

// Driver for routine laguer

int main(void)
{
  const int M=4; // degree of polynomial
  const int MP1=M+1; // no. of coefficients
  const int NTRY=21;
  const double EPS=1.0e-6;
  const complex<double> real1(1.0,0.0), imag1(0.0,1.0);

  // A series of trial values along the complex plane from -1-i to 1+i
  const complex<double> a_d[MP1] =
    {2.0*imag1,0.0,-real1-2.0*imag1,0.0,real1};

  bool newroot;
  int i, its,j,n=0;
  complex<double> x;
  VecComplex a(MP1), y(NTRY);

  for (i=0;i<MP1;i++) a[i] = a_d[i];

  cout << endl << "Roots of polynomial x^4 - (1+2i)*x^2 + 2i" << endl;
  cout << endl << setw(22) << "Root";
  cout << setw(15) << "#iter" << endl << endl;
  cout << fixed << setprecision(6);

  for (i=0;i<NTRY;i++) {
    x = complex<double>((i-10.0)/10.0,(i-10.0)/10.0);
    laguer(a,x,its);
    if (n == 0) {
      n=1;
      y[0]=x;
      cout << setw(5) << n << setw(25) << x;
      cout << setw(6) << its << endl;
    } else {
      newroot=true;
      for (j=0;j<n;j++)
	// Check to see if this root has already been found
	if (abs(x-y[j]) <= EPS*abs(x)) newroot=false;
      if (newroot) {
	y[n++]=x;
	cout << setw(5) << n << setw(25) << x;
	cout << setw(6) << its << endl;
      }
    }
  }
  return 0;
}
