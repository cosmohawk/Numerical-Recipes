#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "../code/nr3.h"
#include "mysort.h"
#include "myran.h"

using namespace std;

// Driver for routine sort3

int main(void) 
{
  int i;
  string txt;
  string amsg, bmsg, cmsg;
  Ran myran(30051986);

  // Prompt for message to scramble
  cout << "Enter message" << endl;
  getline(cin,amsg);
  cout << endl;

  // measure length of string
  const int NLEN = amsg.length(); 
  VecDoub rnds(NLEN), b(NLEN), c(NLEN);

  // generate NLEN random numbers
  for (i=0;i<NLEN;i++) rnds[i] = myran.doub();
 
  // create arrays b and c
  for (i=0;i<NLEN;i++) {
    b[i]=i;
    c[i]=NLEN-1-i;
  }
  
  // print original message
  cout << "original message:" << endl << amsg << endl;
  
  // sort array a while mixing b and c
  sort3(rnds,b,c);
  
  // scramble message according to array b and print.
  bmsg=amsg;
  for (i=0;i<NLEN;i++) bmsg[i]=amsg[int(b[i])];
  cout << endl << "scrambled message:" << endl << bmsg << endl;

  // unscrmable according to array c
  cmsg=amsg;
  for (i=0;i<NLEN;i++) cmsg[int(c[i])] = bmsg[i];
  cout << endl << "mirrired message:" << endl << cmsg << endl;
  
  return 0;
}
