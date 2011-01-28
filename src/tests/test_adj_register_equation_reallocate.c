#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_data_structures.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_test_main.h"

void test_adj_register_equation_reallocate(void)
{
  int cnt, i;
  adj_adjointer adjointer;
  adj_variable var;
  adj_block block;
  adj_equation equation;

  adj_create_adjointer(&adjointer);

  adj_create_block("IdentityOperator", NULL, NULL, 0, &block);
  for (i = 0; i < 2*ADJ_PREALLOC_SIZE + 1; i++)
  {
    adj_create_variable("Velocity", i, 0, ADJ_NORMAL_VARIABLE, &var);
    adj_create_equation(var, 1, &block, &var, &equation); 
    adj_register_equation(&adjointer, equation);
    adj_destroy_equation(&equation);
  }

  adj_equation_count(&adjointer, &cnt);
  adj_test_assert(cnt == 2*ADJ_PREALLOC_SIZE+1, "Wrong equation count");

  adj_destroy_adjointer(&adjointer); 
}
