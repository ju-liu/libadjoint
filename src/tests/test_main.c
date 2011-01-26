#include "adj_setup.h"

int main(void) 
{
  int ierr=0;
  adj_init(ierr);
  TESTNAME();
  adj_finalize(ierr);
  return 0;
}

