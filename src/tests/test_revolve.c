#include <stdlib.h>
#include "libadjoint/revolve_c.h"
#include "libadjoint/adj_test_main.h"

void test_revolve(void)
{
  int snaps;
  struct Revolve *r = NULL;

  snaps = revolve_adjust(r, 10);
  adj_test_assert(snaps == 3);

  r = revolve_create_offline(10, 10);
  revolve_destroy(r);

}
