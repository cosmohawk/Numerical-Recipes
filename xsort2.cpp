#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "mysort.h"
#include "myprint_array.h"

using namespace std;

// Driver for routine sort2

int main(void)
{
  const int NP=100;
  string txt;
  int i;
  VecDoub a(NP), b(NP);
  string fn;

  // take file name as input
  cout << "Please enter data file name." << endl;
  cin >> fn;

  // read in datafile
  ifstream fp(fn);
  if (fp.fail()) {
    string em = "data file " + fn + " not found";
    throw(em.c_str());
  }
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];

  // generate b-array
  for (i=0; i<NP; i++) b[i] = i;

  // sort a and mix b
  sort2(a,b);
  cout << endl << "After sorting a and mixing b, array a is:" << endl;
  cout << fixed << setprecision(2);
  print_array(a,10,7);
  cout << endl << "press return to continue.." << endl;
  cin.get();

  // sort b and mix a
  sort2(b,a);
  cout << endl << "After sorting b and mixing a, array a is:" << endl;
  print_array(b,10,7);
  return 0;
}
  

