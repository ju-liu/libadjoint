#include "adj_adjointer_routines.h"
#include "adj_data_structures.h"
#include "adj_test_tools.h"

void test_adj_forget_equation()
{
  int ierr, cnt;
  adj_adjointer adjointer;
  adj_create_adjointer(&adjointer);
  adj_variable u[5];
  adj_block blocks[10];
  adj_equation equation;

  adj_init_error_codes();

  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &u[0]);
  adj_create_variable("Velocity", 1, 0, ADJ_NORMAL_VARIABLE, &u[1]);
  adj_create_variable("Velocity", 1, 1, ADJ_NORMAL_VARIABLE, &u[2]);
  adj_create_variable("Velocity", 2, 0, ADJ_NORMAL_VARIABLE, &u[3]);
  adj_create_variable("Velocity", 2, 1, ADJ_NORMAL_VARIABLE, &u[4]);

  adj_create_block("IdentityOperator", NULL, NULL, 0, &blocks[0]);
  ierr=adj_create_equation(u[0], 1, &blocks[0], &u[0], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==0, "Should have worked");

  adj_create_nonlinear_block("AdvectionOperator", 1, &u[0], 0.5, NULL, &V);
  adj_create_block("TimesteppingOperator", NULL, NULL, 0, &blocks[1]);
  adj_create_block("BurgersOperator", NULL, NULL, 0, &blocks[2]);
  ierr=adj_create_equation(u[1], 2, &blocks[1], &u[0], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_register_equation(&adjointer, equation);
  adj_chkierr(ierr);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==0, "Should have worked");

  adj_create_block("TimesteppingOperator", NULL, NULL, 0, &blocks[3]);
  adj_create_block("BurgersOperator", NULL, NULL, 0, &blocks[4]);
  ierr=adj_create_equation(u[2], 2, &blocks[3], &u[1], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==0, "Should have worked");

  adj_create_block("TimesteppingOperator", NULL, NULL, 0, &blocks[5]);
  adj_create_block("BurgersOperator", NULL, NULL, 0, &blocks[6]);
  ierr=adj_create_equation(u[3], 2, &blocks[5], &u[2], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_register_equation(&adjointer, equation);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_test_assert(ierr==0, "Should have worked");

  adj_create_block("TimesteppingOperator", NULL, NULL, 0, &blocks[7]);
  adj_create_block("BurgersOperator", NULL, NULL, 0, &blocks[8]);
  ierr=adj_create_equation(u[4], 2, &blocks[7], &u[3], &equation);
  adj_test_assert(ierr==0, "Should have worked");
  adj_chkierr(ierr);
  ierr=adj_register_equation(&adjointer, equation);
  adj_chkierr(ierr);
  adj_test_assert(ierr==0, "Should have worked");
  ierr=adj_destroy_equation(&equation);
  adj_chkierr(ierr);
  adj_test_assert(ierr==0, "Should have worked");



  adj_destroy_adjointer(&adjointer); 
}
