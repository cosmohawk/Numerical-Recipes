#include <iostream>
#include <iomanip>
#include <cmath>
#include "../code/nr3.h"
#include "../code/bessel.h"
#include "mymins.h"

using namespace std;

// Driver for routine golden.

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

int main(void)
{
  const Doub TOL=1.0e-6,EQL=1.0e-3;
  bool newroot;
  int i, j, m, nmin=0;
  Doub ax, bx, cx, fa, fb, fc, xmin;
  VecDoub amin(20);

  Funcd fx;

  cout << "Please select bracketing minimisation method:" << endl;
  cout << setw(6) << "1)" << setw(40) << "Golden ratio" << endl;
  cout << setw(6) << "2)" << setw(40) << "Brent's method" << endl;
  cout << setw(6) <<"3)"<<setw(40)<<"Brent's method (with derivative)"<<endl;

  cin >> m;

  Golden golden;
  Brent brent;
  Dbrent dbrent; 

  cout << "Minima of the function bessj0" << endl;
  cout << setw(10) << "min. #" << setw(9) << "x";
  cout << setw(18) << "bessj0(x)" << setw(13) << "bessj1(x)" << endl;
  cout << fixed << setprecision(6);

  for (i=0;i<100;i++) {

    ax = double(i);
    bx = double(i+1);
    
    // Try new bracketing interval and locate minimum
    switch (m) {
    case 1:
      golden.bracket(ax,bx,fx);
      xmin=golden.minimize(fx);
      break;
    case 2:
      brent.bracket(ax,bx,fx);
      xmin=brent.minimize(fx);
      break;
    case 3:
      dbrent.bracket(ax,bx,fx);
      xmin=dbrent.minimize(fx);
      break;
    }

    // Check if minimum is new.
    // j1 should be 0 at extrema of j0
    if (nmin == 0) {
      amin[nmin++] = xmin;
      cout << setw(7) << nmin << setw(16) << xmin;
      cout << setw(13) << fx(xmin);
      cout << setw(13) << fx.df(xmin);
      cout << endl;
    } else {
      newroot = true;
      for (j=0;j<nmin;j++)
	if (fabs(xmin-amin[j]) <= EQL*xmin) newroot = false;
      if (newroot) {
	amin[nmin++] = xmin;
	cout << setw(7) << nmin << setw(16) << xmin;
	cout << setw(13) << fx(xmin);
	cout << setw(13) << fx.df(xmin);
	cout << endl;
      }
    }
  }
  return  0;
}

    
