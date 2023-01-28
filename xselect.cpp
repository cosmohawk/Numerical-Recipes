#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "mysort.h"
#include "../code/selip.h"

using namespace std;

// Driver for routine selct.

int main(void)
{
  const int NP=100;
  string txt;
  int i,k;
  double q,s;
  VecDoub a(NP),b(NP);
  string fn;

  // Prompt for data and load;
  cout << "Please enter data file name" << endl;
  cin >> fn;
  cin.ignore(); // ignore enter key
  ifstream fp(fn);
  if (fp.fail()) {
    string err_msg = "File not found: " +fn;
    throw(err_msg.c_str());
  }
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];
  
  cout << fixed << setprecision(2);
  for (;;) {
    cout << endl << "Input k (negative to end)" << endl;
    cin >> k;
    if (k <0 || k >= NP) break;
    for (i=0;i<NP;i++) b[i]=a[i];
    // compare two selection algorithms
    s=selip(k,a);
    q=select(k,b);
    cout << "Element at sort index " << k;
    cout << " is " << setw(6) << q << endl;
    cout << "Cross-check from selip routine " << setw(6) << s << endl;
  }
  cout << "Normal completion" << endl;
  return 0;
}
