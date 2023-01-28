#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "../code/nr3.h"
#include "../code/sort.h"

using namespace std;

// Driver for routine hpsel.

int main(void) 
{
  const int NP=100;
  string txt, fn;
  int i,j,h;
  double check;
  VecDoub a(NP);

  // Prompt for data file and load
  cout << "Please enter data file name" << endl;
  cin >> fn;
  cin.ignore();
  ifstream fp(fn);
  if (fp.fail()) {
    string err_msg = "File not found: "+fn+"\n";
    throw(err_msg.c_str());
  }
  getline(fp,txt);
  for (i=0;i<NP;i++) fp >> a[i];

  cout << fixed << setprecision(2);
  for (;;) {
    cout << endl << "Input heap size (or -ve to end)" << endl;
    cin >> h;
    if (h < 0 || h >= NP) break;
    Heapselect hpsel(h);
    cout << "1st" << setw(6) << "2nd" << setw(6) << "3rd" << endl;
    for (i=0;i<10;i++) {
      for (j=0;j<10;j++) {
	hpsel.add(a[10*i+j]);
	if (i+j >0  && (10*i+j)%10 == 0) {
	  cout<<hpsel.report(1)<<setw(6);
	  cout<<hpsel.report(2)<<setw(6);
	  cout<<hpsel.report(3)<< endl;
	}
      }
    }
  }

  return 0;
}
