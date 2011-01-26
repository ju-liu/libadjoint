#include <stdio.h>

void adj_test_assert(unsigned char passed, char *testdesc)
{
    if (passed!=0) 
      printf("  fail: %s\n", testdesc);
    else
      printf("  pass: %s\n", testdesc);
}
