#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "../code/bessel.h"
#include "myroots.h"

/* zbrak takes an interval and then subdivides it into N equal parts. It then 
   identifies any of the subintervals that contain roots. Here we look for 
   roots in J_0 between X1=1 and X2=50. */

// Driver for routine zbrak                

double fx(const double x)
{
  Bessjy bess;
  return bess.j0(x);
}

int main(void)
{
  const int N=100,NBMAX=20;
  const double X1=1.0, X2=50.0;
  int i, nb=NBMAX;
  VecDoub xb1(NBMAX),xb2(NBMAX);

  // Do bracketing
  zbrak(fx,X1,X2,N,xb1,xb2,nb);
  
  // Print results to screen
  cout << endl << "Brackets for the roots of bess.j0: " << endl;
  cout << setw(19) << "lower" << setw(11) << "upper";
  cout << setw(16) << "f(lower)" << setw(11) << "f(upper)" << endl;
  cout << fixed << setprecision(4);
  
  for (i=0;i<nb;i++) {
    cout << " root " << setw(2) << (i+1) << setw(11) << xb1[i];
    cout << setw(11) << xb2[i] << "  " << setw(11) << fx(xb1[i]);
    cout << setw(11) << fx(xb2[i]) << endl;
  }
  
  return 0;
}
