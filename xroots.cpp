#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "../code/bessel.h"
#include "myroots.h"

/* rtbis begins with brackets for the root and then finds the root intself by
   bisection. Here we find all the roots of J_0 between 1 and 50. */

struct Funcd {
  // Structure containing the function and its derivative if known
  Bessjy bess;
  Doub operator() (const Doub x) {
    return bess.j0(x);
  }
  Doub df(const Doub x) {
    return -bess.j1(x);
  }
};

void find_roots(int mthd)
{
  // Number of brackets to make
  const int N=100, NBMAX=20;
  
  // x range in which to look for roots.
  const double X1=1.0, X2=50.0;
  
  // The root itself and the accuracy to which we find it.
  double root, xacc;

  int i, nb;
  VecDoub xb1(NBMAX), xb2(NBMAX);

  // Initiate functor object
  Funcd fx;

  // Do the bracketing step.
  zbrak(fx,X1,X2,N,xb1,xb2,nb);

  // Do root finding and print results to screen.
  cout << endl << "Roots of bess.j0:" << endl;
  cout << setw(20) << "x" << setw(16) << "f(x)" << endl << endl;
  cout << fixed << setprecision(6);

  for (i=0;i<nb;i++) {
    xacc = (1.0e-6)*(xb1[i]*xb2[i])/2.0;
    switch (mthd) {
    case 1:
      root = rtbis(fx,xb1[i],xb2[i],xacc);
      break;
    case 2:
      root = rtflsp(fx,xb1[i],xb2[i],xacc);
      break;
    case 3:
      root = rtsec(fx,xb1[i],xb2[i],xacc);
      break;
    case 4:
      root = zriddr(fx,xb1[i],xb2[i],xacc);
      break;
    case 5:
      root = zbrent(fx,xb1[i],xb2[i],xacc);
      break;
    case 6:
      root = rtnewt(fx,xb1[i],xb2[i],xacc);
      break;
    case 7:
      root = rtsafe(fx,xb1[i],xb2[i],xacc);
      break;
    }
    cout << "root " << setw(3) << (i+1) << setw(15) << root;
    cout << setw(15) << fx(root) << endl;
  }
}

int main(void) 
{

  int mthd;

  // Prompt for root finding method
  cout << "Please choose a root finding method :" << endl;
  cout << setw(6) << "1: " << "bisection" << endl;
  cout << setw(6) << "2: " << "false position" << endl;
  cout << setw(6) << "3: " << "secant method" << endl;
  cout << setw(6) << "4: " << "Ridders' method" << endl;
  cout << setw(6) << "5: " << "Brent's method" << endl;
  cout << setw(6) << "6: " << "Newton-Raphson" << endl;
  cout << setw(6) << "7: " << "Newton-Raphson with bisection" << endl;

  cin >> mthd;

  find_roots(mthd);

}
  

  

