#include <iostream>
#include <iomanip>
#include <cmath>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
//#include "/Users/adam/NumericalRecipies/code/ran.h"
#include "myran.h"

using namespace std;

// Driver for routines Ran, Ranq1, Ranq2
// Numerical Recipes Example Book, Ch. 7, p. 126
// needed to make changes to book version
 
double fnc(const double x1, const double x2, const double x3, const double x4)
{
  return sqrt(x1*x1+x2*x2+x3*x3+x4*x4);
}

// a simple template will work with any type that provides the appropriate members
template <typename T>
void integ(T & func)
{
  const unsigned long twotoj[]={0x1L, 0x2L, 0x4L, 0x8L, 0x10L,
				  0x20L, 0x40L, 0x80L, 0x100L,
				  0x200L, 0x400L, 0x800L, 0x1000L,
				  0x2000L, 0x4000L, 0x8000L};
  const double PI=3.141592653589793238;
  int i,idum=(-1),j,jpower,k;
  double x1,x2,x3,x4;
  VecInt iy(3);
  VecDoub yprob(3);

  // Calculates pi statistically using volume of unit n-sphere
  for (i=0;i<3;i++) iy[i]=0;
  cout << "volume of unit n-sphere, n = 2, 3, 4" << endl;
  cout << "# points    pi    (4/3)*pi    (1/2)*pi^2";
  cout << endl << endl;
  cout << fixed << setprecision(6);
  for (j=0;j<15;j++) {
    for (k=twotoj[j];k>=0;k--) {
      x1=func.doub();
      x2=func.doub();
      x3=func.doub();
      x4=func.doub();
      if (fnc(x1,x2,0.0,0.0) < 1.0) ++iy[0];
      if (fnc(x1,x2,x3,0.0) < 1.0) ++iy[1];
      if (fnc(x1,x2,x3,x4) < 1.0) ++iy[2];
    }
    jpower=twotoj[j+1];
    for (i=0;i<3;i++)
      yprob[i]=double(twotoj[i+2])*iy[i]/jpower;
    for (i=0;i<3;i++)
      yprob[i]=double(twotoj[i+2])*iy[i]/jpower;
    cout << setw(6) << jpower << setw(12) << yprob[0];
    cout << setw(12) << yprob[1] << setw(12) << yprob[2] << endl;
  }
  cout << endl << "actual" << setw(12) << PI;
  cout << setw(12) << 4.0*PI/3 << setw(12) << 0.5*PI*PI << endl;
  cout << "Press enter" << endl;
  cin.get();
}

int main(void)
{
  Ran myran(30051986);
  cout << endl << "Testing Ran: " << endl;
  integ(myran);

  Ranq1 myranq1(30051986);
  cout << endl << "Testing Ranq1: " << endl;
  integ(myranq1);

  Ranq2 myranq2(30051986);
  cout << endl << "Testing Ranq2: " << endl;
  integ(myranq2);

  Ranbyte myranbyte(30051986);
  cout << endl << "Testing Ranbyte: " << endl;
  integ(myranbyte);

  Ranhash myranhash;
  cout << endl << "Testing Ranhash: " << endl;
  cout << "double : " << "\t\t" << myranhash.doub(30051986) << endl;
  cout << "64-bit integer : " << "\t" << myranhash.int64(30051986) << endl;

}
