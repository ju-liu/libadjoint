#include "libadjoint/adj_data_structures.h"
#include "libadjoint/adj_test_tools.h"

void test_adj_variable_str(void)
{
  char buf[255], truebuf[255];
  adj_variable u;

  memset(truebuf, 0, 255 * sizeof(char));
  truebuf[0] = 'V'; truebuf[1] = ':'; truebuf[2] = '0';
  truebuf[3] = ':'; truebuf[4] = '0'; truebuf[5] = ':';
  truebuf[6] = 'F'; truebuf[7] = 'o'; truebuf[8] = 'r';
  truebuf[9] = 'w'; truebuf[10] = 'a'; truebuf[11] = 'r';
  truebuf[12] = 'd';

  adj_create_variable("V", 0, 0, ADJ_NORMAL_VARIABLE, &u);
  adj_variable_str(u, buf, 255);
  adj_test_assert(memcmp(buf, truebuf, 255*sizeof(char)) == 0, "Should be exactly equivalent");
}
