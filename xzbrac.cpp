#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "../code/bessel.h"
#include "myroots.h"

/* zbrac is a root-bracketing routine that works by expanding the range of an
   interval geometrically until it brackets a root. */

// Driver for routine zbrac

double fx(const double x)
{
  Bessjy bess;
  return bess.j0(x);
}

int main(void)
{
  bool success; 
  int i;
  double x1,x2;

  cout << setw(21) << "bracketing values:";
  cout << setw(24) << "function values:" << endl << endl;
  cout << setw(9) << "x1" << setw(12) << "x2";
  cout << setw(21) << "bess.j0(x1)" << setw(13) << "bess.j0(x2)" << endl;
  cout << fixed << setprecision(6);
  
  for (i=0;i<10;i++) {
    x1 = double(i) + 1.0;
    x2 = x1 + 1.0;
    success = zbrac(fx,x1,x2);
    if (success) {
      cout << setw(12) << x1 << setw(12) << x2;
      cout << setw(5) << " " << setw(12) << fx(x1);
      cout << setw(13) << fx(x2) << endl;
    }
  }
  return 0;
}
