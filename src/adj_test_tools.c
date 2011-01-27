#include <stdio.h>
#include "libadjoint/adj_test_tools.h"

void adj_test_assert(unsigned char passed, char *testdesc)
{
  if (!passed) 
    printf("  fail: %s\n", testdesc);
  else
    printf("  pass\n");
}
