#include <string>
#include <iostream>
#include <iomanip>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/code/gamma.h"
#include "myran.h"
//#include "mydeviates.h"
#include "myprimpolytest.h"

using namespace std;

// Implement the test for primitive polynomials.

int main(void){

  // The input polynomial and the result
  Ullong iin;
  int res;
  int i;

  // Initiate test structure
  Primpolytest mytest;

  // Prompt user for input
  cout << "Please enter a 32-bit integer: " << endl;
  cin >> iin;

  // Do test.
  res = mytest.test(iin);

  // Print result.
  if (res == 1) {
    cout << iin << " is a primitive polynomial." << endl;
  } else {
    cout << iin << " is not a primitive polynomial." << endl;
  }
    
}
