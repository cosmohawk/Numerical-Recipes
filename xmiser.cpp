#include <iostream>
#include <iomanip>
#include <cmath>
#include "../code/nr3.h"
#include "myran.h"
#include "../code/ranpt.h"
#include "../code/miser.h"

using namespace std;

/* Miser is a routine for multidimensional Monte Carlo integration with 
   recursive stratified sampling. This driver performs the same integration
   repeatedly, so as to test whether miser's claimed standard deviation in fact
   corresponds to the actual spread of its answers around the known correct
   value. We integrate the same narrow peaked Gaussian as for Vegas. */

// Driver routine for miser

// Define some global variables
int idum; // for ranpt
int ndim; // for func
double xoff;

// The function to integrate
double func(VecDoub_I &pt) {

  int j;
  double ans, sum;

  for (sum=0.0,j=0;j<ndim;j++) sum += (100.0*SQR(pt[j]-xoff));
  ans = (sum < 80.0 ? exp(-sum) : 0.0);
  ans *= pow(5.64189,double(ndim));
  return ans;

}

int main(void) { 

  int j,n,nt,ntries;
  double ave,dith,sumav,sumsd,var;

  // Prompt for random seed and make it negative.
  cout << "IDUM = " << endl;
  cin >> idum;
  if (idum > 0) idum = -idum;

  for (;;) {

    // Prompt for other variables.
    cout << "Enter N NDIM XOFF DITH and NTRIES (N=0 to stop)" << endl;
    cin >> n;
    if (n<=0) break;
    cin >> ndim >> xoff >> dith >> ntries;

    // Initiate variables
    VecDoub regn(2*ndim);
    sumav=sumsd=0.0;

    // Iterations
    for (nt=0;nt<ntries;nt++) {
      
      // Initiate vectors
      for (j=0;j<ndim;j++) {
	regn[j]=0.0;
	regn[ndim+j]=1.0;
      }

      // Run Miser
      miser(func,regn,n,dith,ave,var);
      sumav += SQR(ave-1.0);
      sumsd += sqrt(fabs(var));
    }

    // Calculate st. dev. and compare to predicted value
    sumav = sqrt(sumav/ntries);
    sumsd /= ntries;
    cout << "Fractional error: actual, indicated = ";
    cout << fixed << setprecision(6);
    cout << setw(12) << sumav << setw(13) << sumsd << endl << endl;
  }

  cout << "Normal completion" << endl;
  return 0;
}

