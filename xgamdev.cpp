#include <string>
#include <iostream>
#include <iomanip>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/code/gamma.h"
#include "myran.h"
#include "mydeviates.h"

using namespace std;

/* Driver for gamma distributed random variable.
   Displays a gamma distribution of order ia specified by the user. */

int main(void)
{
  const int N=20, NPTS=10000,ISCAL=200,LLEN=50;
  int i, j, klim, idum=(-13);
  double dd, ia, ibeta;
  VecInt dist(N+1);

  for (;;) {
    // make array of zeros
    for (j=0;j<=N;j++) dist[j]=0;
    
    // Instruct user to provide input, read in this input.
    do {
      cout << endl << "Select order of Gamma distribution (n=1.." << N;
      cout << "), -1 to end" << endl;
      cin >> ia;
    } while (ia > N);
    
    // The instruction "-1" ends the program.
    if (ia < 0) break;

    cout << endl << "Select scaling parameter, beta. " << endl;
    cin >> ibeta;

    // Otherwise draw numbers from a Gamma distribution
    Gammadev gamdev(ia,ibeta,idum);
    for (i=0;i<NPTS;i++) {
      j = int(gamdev.dev());
      if ((j >= 0) && (j <= N)) ++dist[j];
    }
    
    // Print results to screen.
    cout << "Gamma-ditribution deviate, order " << ia;
    cout << " of " << NPTS << " points " << endl << endl;
    cout << setw(6) << "x" << setw(8) << "p(x)";
    cout << setw(10) << "graph:" << endl << endl;
    cout << fixed << setprecision(4);
    for (j=0;j<N;j++) {
      dd = double(dist[j])/NPTS;
      klim=int(ISCAL*dd);
      if (klim > LLEN) klim=LLEN;
      string txt(klim,'*');
      cout << setw(6) << j << setw(8) << dd << " " << txt << endl;
    }
  }
  return 0;
}

