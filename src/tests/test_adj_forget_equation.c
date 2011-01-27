#include "adj_adjointer_routines.h"
#include "adj_data_structures.h"
#include "adj_test_tools.h"

void test_adj_forget_equation()
{
  int ierr, cnt;
  adj_adjointer adjointer;
  adj_create_adjointer(&adjointer);
  adj_variable vars[5];
  adj_block blocks[10];
  adj_equation equation;

  adj_init_error_codes();

  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &vars[0]);
  adj_create_variable("Velocity", 1, 0, ADJ_NORMAL_VARIABLE, &vars[1]);
  adj_create_variable("Velocity", 1, 1, ADJ_NORMAL_VARIABLE, &vars[2]);
  adj_create_variable("Velocity", 2, 0, ADJ_NORMAL_VARIABLE, &vars[3]);
  adj_create_variable("Velocity", 1, 1, ADJ_NORMAL_VARIABLE, &vars[4]);

  adj_create_block("IdentityOperator", NULL, NULL, 0, &blocks[0]);
  ierr=adj_create_equation(vars[0], 1, &blocks[0], &vars[0], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==0, "Should have worked");

  adj_create_block("TimesteppingOperator", NULL, NULL, 0, &blocks[1]);
  adj_create_block("BurgersOperator", NULL, NULL, 0, &blocks[2]);
  ierr=adj_create_equation(vars[1], 2, &blocks[1], &vars[1], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_register_equation(&adjointer, equation);
  adj_chkierr(ierr);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==0, "Should have worked");
/*
  adj_create_block("TimesteppingOperator", NULL, NULL, 0, &blocks[3]);
  adj_create_block("BurgersOperator", NULL, NULL, 0, &blocks[4]);
  ierr=adj_create_equation(vars[2], 2, &blocks[3], &vars[3], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==0, "Should have worked");

  adj_create_block("TimesteppingOperator", NULL, NULL, 0, &blocks[5]);
  adj_create_block("BurgersOperator", NULL, NULL, 0, &blocks[6]);
  ierr=adj_create_equation(vars[3], 2, &blocks[5], &vars[5], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==0, "Should have worked");

  adj_create_block("TimesteppingOperator", NULL, NULL, 0, &blocks[7]);
  adj_create_block("BurgersOperator", NULL, NULL, 0, &blocks[8]);
  ierr=adj_create_equation(vars[4], 2, &blocks[8], &vars[8], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==0, "Should have worked");
*/

  /*

  ierr=adj_create_equation(vars[0], 1, &blocks[0], &vars[0], &equation); 
  adj_test_assert(ierr==ADJ_ERR_OK, "Should work");
  adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr==ADJ_ERR_OK, "Should work");

  adj_equation_count(&adjointer, &cnt);
  adj_test_assert(cnt == 1, "equation count");

  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr!= ADJ_ERR_OK, "Can't register again");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==ADJ_ERR_OK, "Should work");

  adj_equation_count(&adjointer, &cnt);
  adj_test_assert(cnt == 1, "Wrong equation count");

  ierr=adj_create_equation(vars[1], 2, blocks, vars, &equation); 
  adj_test_assert(ierr== ADJ_ERR_OK, "Should work");
  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr == ADJ_ERR_OK, "Register equation");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==ADJ_ERR_OK, "Should work");

  adj_equation_count(&adjointer, &cnt);
  adj_test_assert(cnt == 2, "Wrong equation count");

  adj_destroy_adjointer(&adjointer); 
  */
}
