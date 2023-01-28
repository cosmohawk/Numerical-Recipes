#include <iostream>
#include <iomanip>
#include <cmath>
#include "../code/nr3.h"
#include "../code/bessel.h"
#include "myamoeba.h"

using namespace std;

// Driver for routine amoeba.

Doub func(VecDoub_I &x)
{
  // An exotic test function. Minimum is at (x,y,z) = (0.5,0.6,0.7)
  Bessjy bess;
  return 0.6-bess.j0(SQR(x[0]-0.5)+SQR(x[1]-0.6)+SQR(x[2]-0.7));
}

int main(void)
{
  const int MP=4,NP=3;
  const Doub FTOL=1.0e-10;
  int i,nfunc,j,choice;
  Doub delta = 1;
  VecDoub x(NP),y(MP),pmin(NP);
  MatDoub p(MP,NP);

  for (i=0;i<MP;i++) {
    for (j=0;j<NP;j++) 
      x[j] = p[i][j]=(i == (j+1) ? 1.0 : 0.0);
    y[i] = func(x);
  }

  Amoeba am(FTOL);

  pmin = am.minimize(y,delta,func);

  cout << "Number of function evaluations: " << am.nfunc << endl;
  cout << "Coordinates of minimum: " << endl;
  cout << setprecision(4);
  cout << "(x,y,z) = " <<" ("<<pmin[0]<<","<<pmin[1]<<","<<pmin[2]<<")";
  cout << endl;
  cout << endl;
  cout << "Value of function at minimum:"<<endl;
  cout << "f(x,y,z) = "<<am.fmin<<endl;

}
  
