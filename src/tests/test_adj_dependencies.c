#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_test_main.h"

void test_adj_dependencies(void)
{
  adj_adjointer adjointer;
  adj_nonlinear_block V;
  adj_block B[2], I;
  adj_variable u[3];
  adj_variable u_tmp[2];
  adj_variable_data* data_ptr;
  adj_equation eqn;
  int ierr;
  int adj_equations_u0[3] = {1, 2, 3};
  int adj_equations_u10[3] = {1, 2, 3};
  int adj_equations_u11[2] = {1, 2};

  adj_create_adjointer(&adjointer);

  /* Set up the problem */
  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &u[0]);
  adj_create_variable("Velocity", 1, 0, ADJ_NORMAL_VARIABLE, &u[1]);
  adj_create_variable("Velocity", 1, 1, ADJ_NORMAL_VARIABLE, &u[2]);

  adj_create_block("IdentityOperator", NULL, NULL, &I);
  ierr = adj_create_equation(u[0], 1, &I, &u[0], &eqn);
  adj_test_assert(ierr == 0, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == 0, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == 0, "Should have worked");

  adj_create_nonlinear_block("AdvectionOperator", 1, &u[0], 0.5, NULL, &V);
  adj_create_block("TimesteppingOperator", &V, NULL, &B[0]);
  adj_create_block("BurgersOperator", &V, NULL, &B[1]);
  ierr = adj_create_equation(u[1], 2, B, u, &eqn);
  adj_test_assert(ierr == 0, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == 0, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == 0, "Should have worked");
  adj_destroy_nonlinear_block(&V);
  adj_destroy_block(&B[0]);
  adj_destroy_block(&B[1]);

  adj_create_nonlinear_block("AdvectionOperator", 2, u, 0.5, NULL, &V);
  adj_create_block("TimesteppingOperator", &V, NULL, &B[0]);
  adj_create_block("BurgersOperator", &V, NULL, &B[1]);
  u_tmp[0] = u[0];
  u_tmp[1] = u[2];
  ierr = adj_create_equation(u[2], 2, B, u_tmp, &eqn);
  adj_test_assert(ierr == 0, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == 0, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == 0, "Should have worked");
  adj_destroy_nonlinear_block(&V);
  adj_destroy_block(&B[0]);
  adj_destroy_block(&B[1]);

  /* OK. Now let's check the dependencies */
  ierr = adj_find_variable_data(&(adjointer.varhash), &(u[0]), &data_ptr);
  adj_test_assert(ierr == 0, "Should have worked");
  adj_test_assert(data_ptr->nadjoint_equations == 3, "Should be necessary for all three adjoint equations");
  adj_test_assert(memcmp(data_ptr->adjoint_equations, adj_equations_u0, 3*sizeof(int)) == 0, "Should be {1, 2, 3}");

  ierr = adj_find_variable_data(&(adjointer.varhash), &(u[1]), &data_ptr);
  adj_test_assert(ierr == 0, "Should have worked");
  adj_test_assert(data_ptr->nadjoint_equations == 3, "Should be necessary for all three adjoint equations");
  adj_test_assert(memcmp(data_ptr->adjoint_equations, adj_equations_u10, 3*sizeof(int)) == 0, "Should be {1, 2, 3}");

  ierr = adj_find_variable_data(&(adjointer.varhash), &(u[2]), &data_ptr);
  adj_test_assert(ierr == 0, "Should have worked");
  adj_test_assert(data_ptr->nadjoint_equations == 2, "Should be necessary for the first two adjoint equations");
  adj_test_assert(memcmp(data_ptr->adjoint_equations, adj_equations_u11, 2*sizeof(int)) == 0, "Should be {1, 2}");
}
