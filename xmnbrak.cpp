#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "../code/bessel.h"
#include "../code/mins.h"

using namespace std;

// Driver for routine the various instances of the Bracketmethod base class.

double fx(const double x)
{
  Bessjy bess;
  return bess.j0(x);
}

int main(void)
{
  int i;
  Doub a,b,ax,bx,cx,fa,fb,fc,xmin;
  
  Bracketmethod brack;

  cout << fixed << setprecision(6);
  cout << setw(14) << "a" << setw(13) << "b";
  cout << setw(13) << "c" << endl;
  for (i=0;i<10;i++) {
    a = i*0.5;
    b = (i+1)*0.5;
    
    brack.bracket(a,b,fx);

    ax = brack.ax;
    bx = brack.bx;
    cx = brack.cx;
    fa = brack.fa;
    fb = brack.fb;
    fc = brack.fc;

    cout << setw(3) << "x" << setw(15) << ax;
    cout << setw(13) << bx << setw(13) << cx << endl;
    cout << setw(3) << "f" << setw(15) << fa;
    cout << setw(13) << fb << setw(13) << fc << endl;
  }
}



