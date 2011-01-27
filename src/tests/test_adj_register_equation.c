#include "adj_adjointer_routines.h"
#include "adj_data_structures.h"
#include "adj_test_tools.h"

void test_adj_register_equation()
{
  int ierr, cnt;
  adj_adjointer adjointer;
  adj_create_adjointer(&adjointer);


  adj_variable *vars=(adj_variable*) malloc(2*sizeof(adj_variable));
  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &vars[0]);
  adj_create_variable("Velocity", 1, 0, ADJ_NORMAL_VARIABLE, &vars[1]);
  adj_block *blocks=(adj_block*) malloc(2*sizeof(adj_block));
  adj_create_block("IdendityOperator", NULL, NULL, 0, &blocks[0]);
  adj_create_block("IdendityOperator", NULL, NULL, 0, &blocks[1]);

  adj_equation equation;
  ierr=adj_create_equation(vars[0], 1, &blocks[0], &vars[1], &equation); /* nonsense */
  adj_test_assert(ierr!=ADJ_ERR_OK, "Should not work");

  ierr=adj_create_equation(vars[0], 1, &blocks[0], &vars[0], &equation); 
  adj_test_assert(ierr==ADJ_ERR_OK, "Should work");
  adj_register_equation(&adjointer, equation);
  ierr=adj_destroy_equation(&equation);
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

  adj_destroy_adjointer(adjointer); 
  free(vars);
  free(blocks);
}
