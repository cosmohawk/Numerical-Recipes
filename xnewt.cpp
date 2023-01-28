#include <iostream>
#include <iomanip>
#include <cmath>
#include "../code/nr3.h"
#include "../code/ludcmp.h"
#include "../code/qrdcmp.h"
#include "myroots_multidim.h"

using namespace std;

/* Tests newt of the non linear system,

   f_0 = x_0^2 + x_1^2 -2,
   f_1 = exp(x_0 - 1) + x_1^3 -2,

   which has roots at (x_0,x_1) = (1,1). mnewt, which lacks a globally 
   convergent strategy, fails on this example. */

// User supplied function

VecDoub funcv(VecDoub_I x) 
{
  VecDoub f(2);
  f[0] = SQR(x[0])+SQR(x[1])-2.0;
  f[1] = exp(x[0]-1.0)+x[1]*SQR(x[1])-2.0;
  return f;
}

// Driver for routine newt

int main(void)
{
  const int N=2;
  bool check;
  int i;
  VecDoub x(N),f(N);
  
  // Initial starting guesses
  x[0]=2.0;
  x[1]=0.5;

  // Run newt
  newt(x,check,funcv);
  
  // check results
  f = funcv(x);
  if (check) cout << "Convergence problems." << endl;
  cout << endl << setw(7) << "Index" << setw(8) << "x";
  cout << setw(13) << "f" << endl << endl;
  cout << fixed << setprecision(6);
  for (i=0;i<N;i++)
    cout << setw(13) << i << setw(13) << x[i] << setw(13) << f[i] << endl;
  return 0;
}
