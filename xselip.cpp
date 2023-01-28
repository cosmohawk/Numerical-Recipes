#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "mysort.h"
#include "../code/selip.h"
#include "myprint_array.h"

using namespace std;

// Driver for routine selip.

int main(void)
{
  const int NP=100;
  string txt;
  int i,k;
  double q;
  VecDoub a(NP), b(NP);
  string fn;

  // Prompt for data file and load
  cout << "Please enter data file name" << endl;
  cin >> fn;
  cin.ignore();
  ifstream fp(fn);
  if (fp.fail()) {
    string err_msg = "File not found : "+fn+"\n";
    throw(err_msg.c_str());
  }
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];

  cout << endl << "original array:" << endl;
  cout << fixed << setprecision(2);
  print_array(a,10,7);
  
  for (;;) {
    cout << endl << "Input k (negative to end)" << endl;
    cin >> k;
    if (k < 0 || k >= NP) break;
    q = selip(k,a);
    cout << "Element at sort index " << k;
    cout << " is " << setw(6) << q << endl;
  }

  cout << "Normal completion" << endl;
  return 0;
}
  
