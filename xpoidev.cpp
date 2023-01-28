#include <string>
#include <iostream>
#include <iomanip>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/code/gamma.h"
#include "myran.h"
#include "mydeviates.h"

using namespace std;

/* Driver for routine Poissondev.
   Produces Poisson distribution with mean xm specified by user. */

int main(void)
{
  const int N=20,NPTS=10000,ISCAL=200,LLEN=50;
  int i,j,klim,idum=(-13);
  double dd, lambda;
  VecInt dist(N+1);

  for (;;) {
    // Make array of zeros.
    for (j=0;j<=N;j++) dist[j]=0;
    // Prompt user for input and process it.
    do {
      cout << endl << "Mean of Poisson distribution (0.0<x<" << N << ".0)";
      cout << "-1 to end: " << endl;
      cin >> lambda;
    } while (lambda > N);
    if (lambda < 0.0) break;

    // Use input to initiate a Poisson distribution.
    Poissondev poi(lambda,idum);
    
    // Collect some deviates and bin them.
    for (i=0;i<NPTS;i++) {
      j = int(0.5+poi.dev());
      if ((j >= 0) && (j <= N)) ++dist[j];
    }

    // Print the results, make the graph.
    cout << fixed << setprecision(4);
    cout << "Poisson-distributed deviate, mean " << lambda;
    cout << "of " << NPTS << " points" << endl;
    cout << setw(5) << "x" << setw(9) << "p(x)";
    cout << setw(11) << "graph:" << endl;
    for (j=0;j<=N;j++) {
      dd=double(dist[j])/NPTS;
      klim=int(ISCAL*dd);
      if (klim > LLEN) klim=LLEN;
      string txt(klim,'*');
      cout << setw(6) << j << setw(9) << dd;
      cout << " " << txt << endl;
    }
    cout << endl;
  }
  return 0;
    
}
