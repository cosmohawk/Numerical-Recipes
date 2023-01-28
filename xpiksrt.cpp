#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include "../code/nr3.h"
#include "myprint_array.h"
#include "mysort.h"

using namespace std;

// Driver for routine piksrt

int main(void)
{
  const int NP=100;
  string txt;
  int i;
  VecDoub a(NP);
  ifstream fp("tarray.dat");

  if (fp.fail())
    throw("Data file tarray.dat not found");
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];
  cout << "original array:" << endl;
  cout << fixed << setprecision(2);
  print_array(a,10,7);
  piksrt(a);
  cout << "sorted array:" << endl;
  print_array(a,10,7);
  return 0;

}
