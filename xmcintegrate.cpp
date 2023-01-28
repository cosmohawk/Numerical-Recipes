#include <iostream>
#include <iomanip>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/code/gamma.h"
#include "myran.h"
#include "mymcintegrate.h"

/* Code to do the MC integration of the torus function. */

int main(void) {

  // Specify the box V
  VecDoub xlo(3), xhi(3);
  xlo[0] = 1.;    xhi[0] = 4.;
  xlo[1] = -3.;   xhi[1] = 4.;
  xlo[2] = -1.;   xhi[2] = 1.;

  /* "Simple" MC integration */
  // Create an instance of MCintegrate
  MCintegrate mymc(xlo, xhi, torusfuncs, torusregion, NULL, 10201);

  // Sample and update answers.
  cout << "Doing Monte Carlo integration of torus piece." << endl;
  mymc.step(1000000);
  mymc.calcanswers();

  // Access and print the answer, and its error.
  cout << "Volume of torus piece = " << mymc.ff[0] << " +/- ";
  cout << mymc.fferr[0] << endl;

  // Do integration for different number of points
  int i;
  int npts;
  for (i=1;i<=10;i++) {
    MCintegrate mymc_loop(xlo, xhi, torusfuncs, torusregion, NULL, 10201);
    npts = pow(3,i);
    mymc_loop.step(npts);
    mymc_loop.calcanswers();
    cout << mymc_loop.ff[0] << " ";
    cout << mymc_loop.fferr[0] << endl;
  }


  /* MC integration with a change of variables */
  cout << endl << "Now our density is a strong function of z: " << endl;
  cout << "rho(x,y,z) = exp(5z)" << endl;
  cout << endl << "Doing MCintegration with change of variables." << endl;
  VecDoub slo(3), shi(3);
  slo[0] = 1.0;                shi[0] = 4.;
  slo[1] = -3.;                shi[1] = 4.;
  slo[2] = 0.2*exp(5.*(-1.));  shi[2] = 0.2*exp(5.*(1.));

  MCintegrate mymc2(slo,shi,torusfuncs,torusregion,torusmap,10201);
  mymc2.step(1000000);
  mymc2.calcanswers();

  //  Access and print the answer, and its error.
  cout << "Volume of torus piece = " << mymc2.ff[0] << " +/- ";
  cout << mymc2.fferr[0] << endl;

}
