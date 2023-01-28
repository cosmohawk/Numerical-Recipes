#include "../code/nr3.h"
#include "../code/psplot.h"
#include "myscrsho.h"
#include "../code/bessel.h"

using namespace std;

/* Demonstrates scrsho by plotting the zero order Bessel function J_0 */

// Driver for routine scrsho

double fx(const double x)
{
  Bessjy bess;
  return bess.j0(x);
}

int main(void)
{
  scrsho(fx);
  return 0;
}
