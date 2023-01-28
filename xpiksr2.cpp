#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "myprint_array.h"
#include "mysort.h"

// Driver for routine piksr2.

int main(void)
{
  const int NP=100;
  string txt;
  int i;
  VecDoub a(NP), b(NP);
  string fn;

  cout << "Please enter matrix filename." << endl;
  cin >> fn;

  ifstream fp(fn);

  if (fp.fail()) {
    string error_message = "Data file '"+fn+"' not found";
    throw(error_message.c_str());
  }
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];
  // generate b-array
  for (i=0;i<NP;i++) b[i]=i;
  // sort a and mix b
  piksr2(a,b);
  cout << endl << "After sorting a and mixing b, array a is:" << endl;
  cout << fixed << setprecision(2);
  print_array(a,10,7);
  cout << endl << "... and array b is:" << endl;
  print_array(b,10,7);
  cin.get();
  // sort b and mix a
  piksr2(b,a);
  cout << endl << "After sorting b and mixing a, array a is:" << endl;
  print_array(a,10,7);
  cout << endl << "... and array b is:" << endl;
  print_array(b,10,7);
  return 0;
}
