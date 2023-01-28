#include <iostream>
#include <iomanip>
#include <cmath>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/code/gamma.h"
#include "myran.h"
#include "mydeviates.h"

using namespace std;

// Driver for routine expdev

int main(void)
{
  // Declare variables.
  const int NPTS=10000;
  const double EE=2.7182818284590;
  int i,j,total=0;
  int idum=(1);
  double expect,xx,y;
  VecInt x(20);
  VecDoub trig(21);

  // Expondev class called with beta and integer seed.
  double beta = 1.0;
  Expondev expdev(beta,idum);
  
  // Create histogram bin edges, make empty bins with zeros.
  for (i=0;i<20;i++) {
    trig[i]=i/20.0;
    x[i]=0;
  }
  trig[20]=1.0;

  // Generate deviates and put them into bins
  for (i=0;i<NPTS;i++) {
    y=expdev.dev();
    for (j=0;j<20;j++)
      if ((y<trig[j+1]) && (y>trig[j])) ++x[j];
  }
  for (i=0;i<20;i++) total += x[i];

  // Print the results.
  cout << endl << "exponential distribution with ";
  cout << NPTS << " points:" << endl;
  cout << "    interval    observed    expected"<< endl << endl;
  cout << fixed << setprecision(6);

  // Compare with theoretical distribution.
  for (i=0;i<20;i++) {
    xx=double(x[i])/total;
    expect=exp(-(trig[i]+trig[i+1])/2.0);
    expect *= (0.05*EE/(EE-1));
    cout << setprecision(2);
    cout << setw(8) << trig[i] << setw(6) << trig[i+1];
    cout << setprecision(6);
    cout << setw(12) << xx << setw(12) << expect << endl;
  }
  return 0;
}
