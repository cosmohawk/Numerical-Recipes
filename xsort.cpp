#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "myprint_array.h"
#include "mysort.h"

using namespace std;

// Driver for routine sort.

int main(void)
{
  const int NP=100;
  string txt;
  int i;
  VecDoub a(NP);
  string fn;

  // Prompt user for data & read in file
  cout << "Please enter name of data file" << endl;
  cin >> fn;
  ifstream fp(fn);
  if (fp.fail()) {
    string em = "Data file '"+fn+"' not found.";
    throw(em.c_str());
  }
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];
  fp.close();

  // Print original array
  cout << endl << "original array:" << endl;
  cout << fixed << setprecision(2);
  print_array(a,10,7);

  // Sort array then print sorted array.
  sort(a);
  cout << endl << "sorted array:" << endl;
  print_array(a,10,7);

  return 0;
}
