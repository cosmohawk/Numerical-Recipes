#include <string>
#include <iostream>
#include <iomanip>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/code/gamma.h"
#include "myran.h"
#include "mydeviates.h"

using namespace std;

// Driver for routine gasdev.

int main(void)
{
  // Declare variables.
  const int N=20,NOVER2=N/2,NPTS=10000,ISCAL=400,LLEN=50;
  int i,j,klim,idum=(-13);
  double dd,x;
  VecInt dist(N+1);
  Normaldev gasdev(0.0,1.0,30051986);

  // Initiate dist as an array of zeros.
  for (j=0;j<=N;j++) dist[j]=0;

  // Generate random points and then put them into bins.
  for (i=0;i<NPTS;i++) {
    x=0.25*N*gasdev.dev();
    j=int(x > 0 ? x+0.5 : x-0.5);
    if ((j >= -NOVER2) && (j <= NOVER2)) ++dist[j+NOVER2];
  }

  // Print results.
  cout << "Normally distributed deviate of " << NPTS;
  cout << "points" << endl << endl;
  cout << setw(5) << "x" << setw(11) << "p(x)";
  cout << setw(10) << "graph:" << endl << endl;
  cout << fixed << setprecision(4);

  // Make a plot.
  for (j=0;j<=N;j++) {
    dd = double(dist[j])/NPTS;
    klim=int(ISCAL*dd);
    if (klim > LLEN) klim=LLEN;
    string txt(klim, '*');
    cout << setw(8) << j/(0.25*N) << setw(9) << dd << " ";
    cout << txt << endl;
  }
  return 0;
}

