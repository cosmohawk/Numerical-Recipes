#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include "../code/nr3.h"
#include "myprint_array.h"
#include "mysort.h"

using namespace std;

int main(void)
{
  const int NP=100;
  string txt;
  int i, j, k, l;
  VecInt irank(NP);
  VecDoub a(NP), b(10);
  string fn;

  // Prompt for data file and load
  cout << "Please enter name of data file." << endl;
  cin >> fn;
  cin.ignore(); // ignore enter
  ifstream fp(fn);
  if (fp.fail()) {
    string err_msg = "Data file '"+fn+"' not found.";
    throw(err_msg.c_str());
  }
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];
  
  // index array and rank index
  Indexx indexx(a);
  indexx.rank(irank);
  //Indexx irank(indexx.indx);

  cout << endl << "original array:" << endl;
  cout << fixed << setprecision(2);
  print_array(a,10,7);
  
  cout << endl << "table of ranks:" << endl;
  print_array(irank,10,7);

  cout << "press return to continue";
  cin.ignore();

  cout << endl << "array sorted according to rank table" << endl;
  for (i=0;i<10;i++) {
    for (j=0;j<10;j++) {
      k=10*i+j;
      for (l=0;l<NP;l++)
	if (irank[l] == k) b[j]=a[l];
    }
    print_array(b,10,7);
  }
  return 0;  

}
  
