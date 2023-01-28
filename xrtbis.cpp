#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "../code/bessel.h"
#include "myroots.h"

/* rtbis begins with brackets for the root and then finds the root intself by
   bisection. Here we find all the roots of J_0 between 1 and 50. */

double fx(const double x)
{
  Bessjy bess;
  return bess.j0(x);
}

int main(void)
{
  // Number of brackets to make
  const int N=100, NBMAX=20;
  
  // x range in which to look for roots.
  const double X1=1.0, X2=50.0;
  
  // The root itself and the accuracy to which we find it.
  double root, xacc;

  int i, nb;
  VecDoub xb1(NBMAX), xb2(NBMAX);

  // Do the bracketing step.
  zbrak(fx,X1,X2,N,xb1,xb2,nb);

  // Do root finding and print results to screen.
  cout << endl << "Roots of bess.j0:" << endl;
  cout << setw(20) << "x" << setw(16) << "f(x)" << endl << endl;
  cout << fixed << setprecision(6);

  for (i=0;i<nb;i++) {
    xacc = (1.0e-6)*(xb1[i]*xb2[i])/2.0;
    root = rtbis(fx,xb1[i],xb2[i],xacc);
    cout << "root " << setw(3) << (i+1) << setw(15) << root;
    cout << setw(15) << fx(root) << endl;
  }

  return 0;
}


  

