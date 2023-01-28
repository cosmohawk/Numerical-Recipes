#include <iostream>
#include <iomanip>
#include "/Users/adam/NumericalRecipies/code/nr3.h"
#include "/Users/adam/NumericalRecipies/code/gamma.h"
#include "myran.h"
#include "myhash.h"

using namespace std;

/* This code declares a hash memory for integers and then stores the birth
   years of some people. */

int main(void) {

  // Declare the hashtable.
  Hash<string, Int, Hashfn2> year(1000,1000);

  // Put some entries in the table.
  year[string("Marie Antoinette")] = 1755;
  year[string("Ludwig van Beethoven")] = 1770;
  year[string("Charles Babbage")] = 1791;

  // How old was Marie Antoinette when Charles Babbage was born?
  Int diff = year[string("Charles Babbage")] -year[string("Marie Antoinette")];
  cout << "Marie Antoinette was "<< diff << " years old ";
  cout << "when Charles Babbage was born." << endl;

  // Instead if using C++strings, we can use null terminated C strings as keys.
  Hash<char, Int, Hashfn2> yearc(1000,1000);
  yearc["Marie Antoinette"[0]] = 1755;
  yearc["Ludwig van Beethoven"[0]] = 1770;
  yearc["Charles Babbage"[0]] = 1791;

  Int diff2 = yearc["Charles Babbage"[0]] -yearc["Marie Antoinette"[0]];

  if (diff2 == diff) {
    cout << "C++ strings and null terminated C strings give the same answer.";
    cout << endl;
  }

  // Store names of people into a hash memory indexed by birth year.
  Mhash<Int,string,Hashfn2> person(1000,1000);
  
  // Put some entries in the table.
  person.store(1775, string("Jane Austen"));
  person.store(1791, string("Charles Babbage"));
  person.store(1767, string("Andrew Jackson"));
  person.store(1791, string("James Buchanan"));
  person.store(1767, string("John Quincy Adams"));
  person.store(1770, string("Ludwig van Beethoven"));
  person.store(1791, string("Samuel Morse"));
  person.store(1755, string("Marie Antoinette"));
  
  // Loop over the years and print the people who were born that year.
  // This example only works with C++ strings.
  string str;
  for (Int i=1750;i<1800;i++) {
    if (person.getinit(i)) {
      cout << "\n" << "born in " << i << ":\n";
      while (person.getnext(str)) cout << str.data() << "\n";
    }
  }

}
