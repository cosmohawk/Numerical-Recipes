#include <iostream>
#include <iomanip>
#include <cmath>
#include "../code/nr3.h"
#include "../code/bessel.h"
#include "mymins.h"
#include "mymins_ndim.h"

using namespace std;

// Driver for routine powell

Doub func(VecDoub_I &x)
{
  Bessjy bess;
  return 0.5-bess.j0(SQR(x[0]-1.0)+SQR(x[1]-2.0)+SQR(x[2]-3.0));
}

struct Func {
  Doub operator()(VecDoub_I &x);
};

int main(void)
{
  const int NDIM=3;
  const Doub FTOL=1.0e-6;
  int i,j,iter;
  Doub fret;
  VecDoub p(NDIM);
  p[0]=1.5;p[1]=1.5;p[2]=2.5;
  MatDoub xi(NDIM,NDIM);

  for (i=0;i<NDIM;i++)
    for (j=0;j<NDIM;j++)
      xi[i][j]=(i ==j ? 1.0 : 0.0);
  
  //  Powell<Func> powell(func);
  
  //powell.minimize(xi);

  //  cout << powell.fret << endl;
  //  cout << powell.iter << endl;
	
}
