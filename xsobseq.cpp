#include <iostream>
#include <iomanip>
#include <fstream>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "mysobseq.h"

using namespace std;

/* The routine sobseq returns successive values of a multidimensional quasi-
   random Sobol' sequence. The following example shows how the routine is 
   initiated via a negative argument, and then useed to obtain the first 2^n
   values of an m-dimensional sequence. */

// Driver for routine sobseq

int main(void) {

  // Declare variables.
  int m;
  int i,n1=(-1);
  int n;
  string fname;
  ofstream myfile;

  // Prompt for input
  cout << "Please enter number of dimensions n: " << endl;
  cin >> m;

  cout << "Please enter n, number of points = 2^n: " << endl;
  cin >> n;

  cout << "Please enter name of output file: " << endl;
  cin >> fname;

  // Initiate sobseq.
  VecDoub x(m);
  sobseq(n1,x);

  myfile.open(fname);

  myfile << fixed << setprecision(5);
  // cout << fixed << setprecision(5);
  for (i=0;i<pow(2,n);i++) {
    sobseq(m,x);
    myfile << setw(11) << x[0] << setw(11) << x[1];
    myfile << setw(11) << x[2] << setw(11) << (i+1) << endl;
    //cout << setw(11) << x[0] << setw(11) << x[1];
    //cout << setw(11) << x[2] << setw(11) << (i+1) << endl;
  }
  myfile.close();
  return 0;
}
