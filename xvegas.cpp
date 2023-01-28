#include <iostream>
#include <iomanip>
#include <cmath>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/mycode/myran.h"
#include "/Users/adam/NumericalRecipies/code/rebin.h"
#include "/Users/adam/NumericalRecipies/code/vegas.h"

using namespace std;

// Driver for routine vegas.

int idum;        // for ranno
int ndim;        // for fxn
double xoff;

double fxn(VecDoub_I &pt, const double wgt) {
  /* A narrow peaked n-dimensional Gaussian function. */

  int j;
  double ans, sum;

  for (sum=0.0,j=0;j<ndim;j++) sum += (100.0*SQR(pt[j]-xoff));
  ans = (sum < 80.0 ? exp(-sum) : 0.0);
  ans *= pow(5.64189,double(ndim));

  return ans;
}

int main(void) {

  // Declare variables.
  int init, itmax, j, ncall, nprn;
  double avgi, chi2a, sd, xoff;

  // Prompt for random seed.
  cout << "IDUM = " << endl;
  cin >> idum;

  // Random seed should be negative.
  if (idum > 0) idum = -idum;
  
  // Prompt for other parameters.
  cout << fixed << setprecision(6);
  for (;;) {
    cout << "Enter NDIM XOFF NCALL ITMAX NPRN (NDIM=0 to stop)" << endl;
    cin >> ndim; // The number of dimensions.
    if (ndim <= 0) break;
    cin >> xoff; // A scalar offset, used to construct vector x0
    cin >> ncall; // The maximum numer of function evaluations per iteration
    cin >> itmax; // The number of adaptive iterations.
    cin >> nprn; // Parameter determining verbosity of vegas.

    VecDoub regn(2*ndim);
    avgi=sd=chi2a=0.0;

    for (j=0;j<ndim;j++) {
      regn[j] = 0.0;
      regn[j+ndim] = 1.0;
    }

    // Run Vegas
    init = -1;
    vegas(regn,fxn,init,ncall,itmax,nprn,avgi,sd,chi2a);
    cout << "Number of iterations performed: " << itmax << endl;
    cout << "Integral, Stndard Dev., Chi-sq. = ";
    cout << setw(13) << avgi << setw(13) << sd;
    cout << setw(13) << chi2a << endl;

    init = 1;
    cout << "Additional iterations performed: " << itmax << endl;
    cout << "Integral, Standard Dev., Chi-sq. = ";
    cout << setw(13) << avgi << setw(13) << sd;
    cout << setw(13) << chi2a << endl << endl;
  }

  cout << "Normal completion" << endl;
  return 0;
}
