#include <iostream>
#include <iomanip>
#include <complex>
#include "../code/nr3.h"
#include "myroots_poly.h"

using namespace std;

/* zroots is a driver for laguer. This routine exercises zroots by first 
   finding the roots to the polynomial, then corrupting them, then polishing
   them. 

   See also xlaguer. */

// Driver for routine zroots.

int main(void)
{
  const int M=4, MP1=M+1;
  const complex<double> real1(1.0,0.0), imag1(0.0,1.0);
  const complex<double> a_d[MP1] = 
    {2.0*imag1,0.0,-real1-2.0*imag1,0.0,real1};
  int i;
  bool polish;
  VecComplex a(MP1), roots(M);

  for (i=0;i<MP1;i++) a[i] = a_d[i];

  // Find unpolished roots of polynomial
  cout << endl <<"Roots of the complex polynomial x^4 - (1+2i)*x^2 + 2i"<<endl;
  polish=false;
  zroots(a,roots,polish);
  cout << endl << "Unpolished roots:" << endl;
  cout << setw(14) << "root #" << setw(14) << "root:" << endl << endl;
  cout << fixed << setprecision(6);
  for (i=0;i<M;i++) {
    cout << setw(11) << noshowpos << i;
    cout << setw(25) << showpos << roots[i] << endl;
  }
  
  // Corrupt roots
  cout << endl << "Corrupted roots:" << endl;
  for (i=0;i<M;i++) 
    roots[i] = (double(1.0)+double(0.01)*(i+1))*roots[i];
  cout << setw(14) << "root #" << setw(14) << "root :" << endl << endl;
  for (i=0;i<M;i++) {
    cout << setw(11) << noshowpos << i;
    cout << setw(25) << showpos << roots[i] << endl;
  }
  
  // Repolish roots
  polish=true;
  zroots(a,roots,polish);
  cout << endl << "Polish roots:" << endl;
  cout << setw(14) << "root #" << setw(14) << "root:" << endl << endl;
  for (i=0;i<M;i++) {
    cout << setw(11) << noshowpos << i;
    cout << setw(25) << showpos << roots[i] << endl;
  }
  return 0;
  
}
