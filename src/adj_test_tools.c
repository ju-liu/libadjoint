#include <stdio.h>
#include "libadjoint/adj_test_tools.h"

void adj_test_assert(int passed, char *testdesc)
{
  if (!passed) 
    printf("  fail: %s\n", testdesc);
  else
    printf("  pass\n");
}

int adj_sizeof_adjointer()
{
  return sizeof(adj_adjointer);
}
