#include <string>
#include <iostream>
#include <iomanip>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/code/gamma.h"
#include "/Users/adam/NumericalRecipies/code/cholesky.h"
#include "myran.h"
#include "mydeviates.h"
#include "mymultinormaldev.h"

using namespace std;

int main(void)
{
  const int NPTS = 10000;
  int i;

  // Declare covariance matrix and mean vector
  MatDoub covmat(2,2);
  VecDoub meanvar(2);

  // some values to put in the matrix
  double cov_aa = 1.2;
  double cov_ab = 1.0;
  double cov_bb = 0.9;

  // Random seed
  int seed = 19021994;

  meanvar[0] = -1.0;
  meanvar[1] = -1.0;

  covmat[0][0] = cov_aa;
  covmat[0][1] = cov_ab;
  covmat[1][1] = cov_bb;

  Multinormaldev mndev(seed, meanvar, covmat);

  for (i=0;i<NPTS;i++) {
    cout << mndev.dev()[0] << " " << mndev.dev()[1] << endl;
  }

  return 0;
}
