#include "adj_setup.h"

int main(void) 
{
  int ierr;
  adj_init(ierr);
  TESTNAME();
  adj_finalise(ierr);
}

