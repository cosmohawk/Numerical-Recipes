#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "../code/sort.h"
#include "myprint_array.h"

using namespace std;

// Driver for routine indexx

int main(void)
{
  const int NP=100;
  string txt;
  int i,j;
  //VecInt indx(NP);
  VecDoub a(NP);
  string fn;

  // Prompt for input data and load
  cout << "Please enter data file name" << endl;
  cin >> fn;
  ifstream fp(fn);
  if (fp.fail()) {
    string error_message = "Data file '"+fn+"' not found";
    throw(error_message.c_str());
  }
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];
  
  // Index array and print results.
  Indexx indexx(a);
  cout << endl << "original array:" << endl;
  cout << fixed << setprecision(2);
  print_array(a,10,7);
  cout << "sorted array:" << endl;
  for (i=0;i<10;i++) {
    for (j=0;j<10;j++) cout << setw(7) << a[indexx.indx[10*i+j]];
    cout << endl;
  }
  return 0;
}


  
