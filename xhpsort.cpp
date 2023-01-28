#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "mysort.h"
#include "myprint_array.h"

using namespace std;

// Driver for routine hpsort

int main(void)
{
  const int NP=100;
  string txt;
  int i;
  VecDoub a(NP);
  string fn;

  // Prompt for data and load
  cout << "Please enter data filename" << endl;
  cin >> fn;
  ifstream fp(fn);
  if (fp.fail()) {
    string em = "Data file " + fn + " not found";
    throw(em.c_str());
  }
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];
  cout << endl << "original array:" << endl;
  cout << fixed << setprecision(2);
  print_array(a,10,7);

  // Do heap sort
  hpsort(a);
  cout << endl << "sorted array:" << endl;
  print_array(a,10,7);
  return 0;
}


