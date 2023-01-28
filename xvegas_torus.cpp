#include <iostream>
#include <iomanip>
#include <cmath>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/mycode/myran.h"
#include "/Users/adam/NumericalRecipies/code/rebin.h"
#include "/Users/adam/NumericalRecipies/code/vegas.h"

using namespace std;

/* Perform the density integral of the torus example using Vegas Monte Carlo 
   integration. 

   The additiontional argument wgt in Doub allows us to perform an additional
   integration of a similar function at the same time. See book for
   explaination. */

Doub torusfunc(const VecDoub &x, const Doub wgt) {
  Doub den = exp(5.*x[2]);
  if (SQR(x[2])+SQR(sqrt(SQR(x[0]) + SQR(x[1]))-3.) <= 1.) return den;
  else return 0;
}

int main(void) {

  Doub tgral, sd, chi2a;
  VecDoub regn(6);
  regn[0] = 1.;
  regn[1] = -3.;
  regn[2] = -1.; 
  regn[3] = 4.;
  regn[4] = 4.;
  regn[5] = 1.;

  int init, ncall, itmax, nprn;
  
  vegas(regn,torusfunc,init=0,ncall=10000,itmax=10,nprn=0,tgral,sd,chi2a);
  //cout<< "tgral = "<< tgral << " sd = " << sd << " chi2a = "<< chi2a << endl;
  vegas(regn,torusfunc,1,900000,1,0,tgral,sd,chi2a);
  //cout<< "tgral = "<< tgral <<" sd = " << sd << " chi2a = "<<chi2a << endl;

}

