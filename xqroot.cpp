#include <iostream>
#include <iomanip>
#include <cmath>
#include "../code/nr3.h"
#include "../code/poly.h"
#include "../code/qroot.h"
//#include "myqroot.h"

using namespace std;

/*qroot is used for finding quadratic factors of polynomials with real 
  coefficients. The example polynomial is:

  P(x) = x^6 -6x^5 + 16x^4 -24x^3 +25x^2 - 18x + 10

  The answer should be:

  P(x) = (x^2 + 1)(x^2 -4x + 5)(x^2 -2x +2).
*/

// Driver for routine qroot

int main(void)
{
  const int N=6; // degree of polynomial
  const int NTRY=10; // Number of iterations
  const double EPS=1.0e-6,TINY=1.0e-5;
  const double p_d[N+1]={10.0,-18.0,25.0,-24.0,16.0,-6.0,1.0};
  bool newroot;
  int i,j,nroot=0;
  VecDoub p(N+1),b(NTRY),c(NTRY);

  for (i=0;i<N+1;i++) p[i]=p_d[i];

  cout << endl << "P(x)=x^6 - 6x^5 + 16x^5 - 24x^3 + 25x^2 - 18x + 10" << endl;
  cout << "Quadratic factors x^2 + bx + c" << endl << endl;
  cout << setw(6) << "factor" << setw(11) << "b";
  cout << setw(13) << "c" << endl << endl;
  cout << fixed << setprecision(6);

  // Try values of the form x^2 + b + c
  // Search along a grid
  for (i=0;i<NTRY;i++) {
    c[i] = 0.5*(i+1);
    b[i] = -0.5*(i+1);
    // Run qroot
    qroot(p,b[i],c[i],EPS);
    if (nroot == 0) {
      cout << setw(4) << nroot << setw(16) << b[i];
      cout << setw(13) << c[i] << endl;
      nroot=1;
    } else {
      // Check to see if root has already been found
      newroot=true;
      for (j=0;j<nroot;j++) 
	if ((fabs(b[i]-b[j]) < TINY) && (fabs(c[i]-c[j]) < TINY))
	  newroot=false;
      if (newroot) {
	cout << setw(4) << nroot << setw(16) << b[i];
	cout << setw(13) << c[i] << endl;
	++nroot;
      }
    }
  }
  return 0;
}
    
