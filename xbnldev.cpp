#include <string>
#include <iostream>
#include <iomanip>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/code/gamma.h"
#include "myran.h"
#include "mydeviates.h"

using namespace std;

// Driver for routine Binomialdev.

int main(void)
{
  // Declare variables.
  const int N=20, NPTS=1000, ISCAL=200, NN=100, LLEN=50;
  int i,j,klim,idum=(-133);
  double pp,xm,dd;
  VecInt dist(N+1);

  for (;;) {
    // Create array of zeros
    for (j=0;j<=N;j++) dist[j]=0;

    // Prompt user for input
    do {
      cout << "Mean of binomial distribution (0.0 to ";
      cout << N << ".0)" << " - Negative to end: " << endl;
      cin >> xm;
      cout << endl;
    } while (xm > N);

    // A negative number ends the program
    if (xm < 0.0) break;

    // Declare binomial distribution
    pp=xm/NN;   
    Binomialdev bnl(NN,pp,idum);

    // Draw values from the distribution.
    for (i=0;i<NPTS;i++) {
      j=bnl.dev();
      if (j >= 0 && j <= N) ++dist[j];
    }

    // Print results and make plot of distribution.
    cout << "Binomial-distributed deviate, mean " << xm << " of ";
    cout << NPTS << " points" << endl;
    cout << setw(4) << "x" << setw(9) << "p(x)";
    cout << setw(11) << "graph:" << endl << endl;
    cout << fixed << setprecision(3);
    for (j=0;j<N;j++) {
      dd = double(dist[j])/NPTS;
      klim=int(ISCAL*dd+1);
      if (klim > LLEN) klim=LLEN;
      string txt(klim,'*');
      cout << setw(4) << j << setw(9) << dd << " " << txt << endl;
    }
    cout << endl;
  } 
  return 0;

}
