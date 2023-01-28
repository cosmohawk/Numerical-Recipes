#include <iostream>
#include <iomanip>
#include <cmath>
#include "../code/nr3.h"
#include "../code/ludcmp.h"
#include "../code/mnewt.h"

using namespace std;

/* mnewt looks for roots of multiple nonlinear equations. The user-supplied 
   function usrfun returns the Jacobian matrix fjac of partial derivatives of
   the function. In this example the equations are:

   0 = -x_0^2 -x_1^2 -x_2^2 + x_3
   0 = x_0^2 + x_1^2 +x_2^2 + x_3^2 - 1
   0 = x_0 - x_1
   0 = x_1 - x_2

*/

// User defined function containing polynomials and Jacobian

void usrfun(VecDoub_I &x, VecDoub_O &fvec, MatDoub_O &fjac)
{
  int i;
  int n=x.size();

  fvec[0] = -SQR(x[0]) - SQR(x[1]) - SQR(x[2]) + x[3]; 
  fvec[1] = SQR(x[0]) + SQR(x[1]) + SQR(x[2]) + SQR(x[3]) -1;
  fvec[2] = x[0] - x[1];
  fvec[3] = x[1] - x[2];

  fjac[0][0] = -2.0*x[0]; 
  fjac[0][1] = -2.0*x[1]; 
  fjac[0][2] = -2.0*x[2];
  fjac[0][3] = -1;
  for (i=0;i<n;i++) fjac[1][i] = 2.0*x[i];
  fjac[2][0] = 1.0;
  fjac[2][1] = -1.0;
  fjac[2][2] = 0.0;
  fjac[2][3] = 0.0;
  fjac[3][0] = 0.0;
  fjac[3][1] = 1.0;
  fjac[3][2] = -1.0;
  fjac[3][3] = 0.0;

}  

// Driver for routine mnewt

int main(void)
{
  const int NTRIAL=5,N=4; // Number of trials and degree of polynomials
  const double TOLX=1.0e-6,TOLF=1.0e-6;
  int kk, k, i, j;
  double xx;
  VecDoub fvec(N),x(N);
  MatDoub fjac(N,N);
  
  // Search to find a set of variables that solve equations
  cout << fixed << setprecision(6);
  for (kk=0;kk<2;kk++) {
    for (k=0;k<3;k++) {
      xx = 0.2001*k*(2*kk-1);
      cout << "Starting vector number " << (k+1) << endl << endl;
      for (i=0;i<4;i++) {
	x[i]=xx+0.2*(i+1);
	cout << setw(7) << "x[" << i << "] = ";
	cout << setw(12) << x[i] << endl;
      }
      cout << endl;
      for (j=0;j<NTRIAL;j++) {
	mnewt(1,x,TOLX,TOLF);
	usrfun(x,fvec,fjac);
	cout << setw(5) << "i" << setw(12) << "x[i]";
	cout << setw(14) << "f" << endl;
	for (i=0;i<N;i++) {
	  cout << setw(5) << i << setw(14) << x[i];
	  cout << setw(15) << fvec[i] << endl;
	}
	cout << endl << "press RETURN to continue..." << endl;
	cin.get();
      }
    }
  }
  return 0;
}





