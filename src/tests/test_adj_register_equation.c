#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_data_structures.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_test_main.h"

void test_adj_register_equation(void)
{
  int ierr, cnt;
  adj_adjointer adjointer;
  adj_create_adjointer(&adjointer);
  adj_variable vars[2];
  adj_block blocks[2];
  adj_equation equation;

  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &vars[0]);
  adj_create_variable("Velocity", 1, 0, ADJ_NORMAL_VARIABLE, &vars[1]);
  adj_create_block("IdentityOperator", NULL, NULL, 1.0, &blocks[0]);
  adj_create_block("IdentityOperator", NULL, NULL, 1.0, &blocks[1]);

  ierr=adj_create_equation(vars[0], 1, &blocks[0], &vars[1], &equation); /* nonsense */
  adj_test_assert(ierr!=ADJ_OK, "Should not work");

  ierr=adj_create_equation(vars[0], 1, &blocks[0], &vars[0], &equation); 
  adj_test_assert(ierr==ADJ_OK, "Should work");
  adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr==ADJ_OK, "Should work");

  adj_equation_count(&adjointer, &cnt);
  adj_test_assert(cnt == 1, "equation count");

  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr!= ADJ_OK, "Can't register again");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==ADJ_OK, "Should work");

  adj_equation_count(&adjointer, &cnt);
  adj_test_assert(cnt == 1, "Wrong equation count");

  ierr=adj_create_equation(vars[1], 2, blocks, vars, &equation); 
  adj_test_assert(ierr== ADJ_OK, "Should work");
  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr == ADJ_OK, "Register equation");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==ADJ_OK, "Should work");

  adj_equation_count(&adjointer, &cnt);
  adj_test_assert(cnt == 2, "Wrong equation count");

  adj_destroy_adjointer(&adjointer); 
}
