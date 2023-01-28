#include <iostream>
#include <iomanip>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/mycode/myhashall.h"

using namespace std;

// Driver for routine myhashall.

int main(void) {
  
  int i;
  VecUint idum(4);
  idum[0] = 1234;
  idum[1] = 8763;
  idum[2] = 9467;
  idum[3] = 7793;

  cout << fixed << setprecision(6);
  hashall(idum);
  cout << endl << "hashall gets values: " << endl;
  for (i=0;i<4;i++) cout << idum[i] << setw(15);

  cout << endl;

  return 0;
}
