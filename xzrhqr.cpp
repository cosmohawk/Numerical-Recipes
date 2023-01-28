#include <iostream>
#include <iomanip>
#include <complex>
#include "../code/nr3.h"
#include "../code/eigen_unsym.h"
#include "myzrhqr.h"
using namespace std;

/* zrhqr is an alternative to laguer for polynomials with real coefficients. */

// Driver for routine zrhqr.

int main(void)
{
  const int M=4; // degree of polynomial.
  const int MP1=M+1; // no. of polynomial coefficients.
  const double a_d[MP1]={-1.0,0.0,0.0,0.0,1.0};
  int i;
  VecDoub a(MP1);
  VecComplex rt(M);

  // Convert a_d to a
  for (i=0;i<MP1;i++) a[i] = a_d[i];

  cout << endl << "Roots of polynomial x^4-1" << endl;
  cout << endl << "    #" << setw(17) << "Root" << endl << endl;
  // Call root finder
  zrhqr(a,rt);
  cout << fixed << setprecision(6);
  for (i=0;i<M;i++) {
    cout << setw(5) << noshowpos << i;
    cout << setw(25) << showpos << rt[i] << endl;
  }
  return 0;
}

  
