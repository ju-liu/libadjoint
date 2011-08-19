#include <stdlib.h>
#include <stdio.h>
#include "libadjoint/revolve_c.h"
#include "libadjoint/adj_test_main.h"
#include "libadjoint/adj_test_tools.h"

void test_revolve(void)
{
  int snaps;
  CACTION action;
  CRevolve r;

  snaps = revolve_adjust(r, 10);
  adj_test_assert(snaps == 3, "Snaps should be 3");

  r = revolve_create_offline(10, 10);
 
  action = revolve(r);
  printf("Action: %i\n", action);
  printf("Action: %s\n", revolve_caction_string(action));

  revolve_destroy(r);

}
