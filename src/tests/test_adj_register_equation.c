#include "adj_adjointer_routines.h"
#include "adj_data_structures.h"
#include "adj_test_tools.h"

void test_adj_register_equation()
{
  int ierr, cnt;
  adj_adjointer adjointer;
  adj_create_adjointer(&adjointer);

  adj_variable a;
  adj_variable b;
  adj_variable c;
  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &a);
  adj_create_variable("Velocity", 1, 0, ADJ_NORMAL_VARIABLE, &b);
  adj_create_variable("Velocity", 2, 0, ADJ_NORMAL_VARIABLE, &c);
  adj_block I;
  adj_create_block("IdendityOperator", NULL, NULL, 0, &I);

  adj_equation equation;
  ierr=adj_create_equation(a, 1, &I, &b, &equation); /* nonsense */
  adj_test_assert(ierr!=0, "Should not work");

  ierr=adj_create_equation(a, 1, &I, &a, &equation); /* should work */
  adj_test_assert(ierr==0, "Should work");

  adj_register_equation(&adjointer, equation);
  adj_equation_count(&adjointer, &cnt);
  adj_test_assert(cnt == 1, "equation count");

  ierr=adj_register_equation(&adjointer, equation); /* can't register again */
  adj_test_assert(ierr!=0, "Should not work");

  adj_equation_count(&adjointer, &cnt);
  printf("%i", cnt);
  adj_test_assert(cnt == 1, "equation count");

  /* adj_destroy_adjointer(adjointer); */
}
