#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "mysort.h"
#include "myprint_array.h"

using namespace std;

// Driver for routine shell.

int main(void)
{
  const int NP=100;
  string txt;
  int i;
  VecDoub a(NP);
  string fn;

  // Prompt user for filename
  cout << "Please enter name of datafile." << endl;
  cin >> fn;
  ifstream fp(fn);
  // Check that file exists.
  if (fp.fail()) {
    string em = "Data file '"+fn+"' not found.";
    throw(em.c_str());
  }
  // Read in file;
  getline(fp,txt);
  for (i=0; i<NP; i++) fp >> a[i];
  fp.close();

  // Print original array
  cout << endl << "Original array: " << endl;
  cout << fixed << setprecision(2);
  print_array(a,10,7);
  
  // Do shell sort and print sorted array
  shell(a);
  cout << endl << "Sorted array: " << endl;
  print_array(a,10,7);
  
  return 0;
}
