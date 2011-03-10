#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_test_main.h"

void test_adj_iteration_count(void)
{
  adj_adjointer adjointer;
  adj_nonlinear_block V;
  adj_block B[2], I;
  adj_variable u[3];
  adj_variable u_tmp[2];
  adj_equation eqn;
  int ierr, count;

  adj_create_adjointer(&adjointer);

  /* Set up the problem */
  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &u[0]);
  adj_create_variable("Velocity", 1, 0, ADJ_NORMAL_VARIABLE, &u[1]);
  adj_create_variable("Velocity", 1, 1, ADJ_NORMAL_VARIABLE, &u[2]);

  adj_create_block("IdentityOperator", NULL, NULL, &I);
  ierr = adj_create_equation(u[0], 1, &I, &u[0], &eqn);
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");

  adj_create_nonlinear_block("AdvectionOperator", 1, &u[0], NULL, &V);
  adj_nonlinear_block_set_coefficient(&V, 0.5);
  adj_create_block("TimesteppingOperator", &V, NULL, &B[0]);
  adj_create_block("BurgersOperator", &V, NULL, &B[1]);
  ierr = adj_create_equation(u[1], 2, B, u, &eqn);
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  adj_destroy_nonlinear_block(&V);
  adj_destroy_block(&B[0]);
  adj_destroy_block(&B[1]);

  adj_create_nonlinear_block("AdvectionOperator", 2, u, NULL, &V);
  adj_nonlinear_block_set_coefficient(&V, 0.5);
  adj_create_block("TimesteppingOperator", &V, NULL, &B[0]);
  adj_create_block("BurgersOperator", &V, NULL, &B[1]);
  u_tmp[0] = u[0];
  u_tmp[1] = u[2];
  ierr = adj_create_equation(u[2], 2, B, u_tmp, &eqn);
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  adj_destroy_nonlinear_block(&V);
  adj_destroy_block(&B[0]);
  adj_destroy_block(&B[1]);

  /* Check that we have indeed two timesteps in the adjointer */
   adj_timestep_count(&adjointer, &count);
   adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
   adj_test_assert(count == 2, "We expect two registered timesteps");

  /* Let's see how many itertion the variables have */
  adj_create_variable("Velocity", 0, -1, ADJ_NORMAL_VARIABLE, &u_tmp[0]);
  ierr = adj_iteration_count(&adjointer, u_tmp[0], &count); 
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  adj_test_assert(count == 1, "Velocity at timestep 0 should have 1 iteration");
  
  adj_create_variable("Velocity", 1, -1, ADJ_NORMAL_VARIABLE, &u_tmp[0]);
  ierr = adj_iteration_count(&adjointer, u_tmp[0], &count); 
  adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked");
  adj_test_assert(count == 2, "Velocity at timestep 0 should have 2 iteration");

  adj_create_variable("Velocity", 2, -1, ADJ_NORMAL_VARIABLE, &u_tmp[0]);
  ierr = adj_iteration_count(&adjointer, u_tmp[0], &count); 
  adj_test_assert(ierr != ADJ_ERR_OK, "Should have not worked");
  adj_test_assert(count == 0, "Velocity at timestep 2 should have 0 iteration");
}
