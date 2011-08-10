#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_test_main.h"

#ifndef HAVE_PETSC
void test_adj_forget_adjoint_equation(void)
{
  adj_test_assert(1 == 1, "Don't have PETSc so can't run this test.");
}
#else
#include "libadjoint/adj_petsc_data_structures.h"

int* has_values(adj_adjointer *adjointer, int nb_vars, adj_variable *vars);

void test_adj_forget_adjoint_equation_fnl(void)
{
  adj_adjointer adjointer;
  adj_nonlinear_block V;
  adj_block B[2], I;
  adj_variable u[5];
  adj_variable u_tmp[2];
  adj_equation eqn;
  int ierr;
  int before_2[5] = {1, 1, 1, 1, 1};
  int after_2[5] = {1, 1, 1, 1, 0};
  int after_1[5] = {1, 1, 1, 0, 0};
  int after_0[5] = {0, 0, 0, 0, 0};
  Vec vec;
  adj_vector value;
  int i, dim=2;
  int* has_val;
  adj_storage_data storage;


  adj_create_adjointer(&adjointer);
  adj_set_petsc_data_callbacks(&adjointer);

  /* Set up the problem */
  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &u[0]);
  adj_create_variable("Velocity", 1, 0, ADJ_NORMAL_VARIABLE, &u[1]);
  adj_create_variable("Velocity", 1, 1, ADJ_NORMAL_VARIABLE, &u[2]);
  adj_create_variable("Velocity", 2, 0, ADJ_NORMAL_VARIABLE, &u[3]);
  adj_create_variable("Velocity", 2, 1, ADJ_NORMAL_VARIABLE, &u[4]);

  adj_create_block("IdentityOperator", NULL, NULL, &I);
  ierr = adj_create_equation(u[0], 1, &I, &u[0], &eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  adj_create_nonlinear_block("AdvectionOperator", 1, &u[0], NULL, &V);
  adj_nonlinear_block_set_coefficient(&V, 0.5);
  adj_create_block("TimesteppingOperator", &V, NULL, &B[0]);
  adj_create_block("BurgersOperator", &V, NULL, &B[1]);
  ierr = adj_create_equation(u[1], 2, B, u, &eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
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
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  adj_destroy_nonlinear_block(&V);
  adj_destroy_block(&B[0]);
  adj_destroy_block(&B[1]);

  adj_create_nonlinear_block("AdvectionOperator", 1, &u[2], NULL, &V);
  adj_nonlinear_block_set_coefficient(&V, 0.5);
  adj_create_block("TimesteppingOperator", &V, NULL, &B[0]);
  adj_create_block("BurgersOperator", &V, NULL, &B[1]);
  ierr = adj_create_equation(u[3], 2, B, &u[2], &eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  adj_destroy_nonlinear_block(&V);
  adj_destroy_block(&B[0]);
  adj_destroy_block(&B[1]);

  adj_create_nonlinear_block("AdvectionOperator", 2, &u[2], NULL, &V);
  adj_nonlinear_block_set_coefficient(&V, 0.5);
  adj_create_block("TimesteppingOperator", &V, NULL, &B[0]);
  adj_create_block("BurgersOperator", &V, NULL, &B[1]);
  u_tmp[0] = u[2];
  u_tmp[1] = u[4];
  ierr = adj_create_equation(u[4], 2, B, u_tmp, &eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  adj_destroy_nonlinear_block(&V);
  adj_destroy_block(&B[0]);
  adj_destroy_block(&B[1]);

  /* Now let's record things, and then forget. */
  VecCreateSeq(PETSC_COMM_SELF, dim, &vec);
  VecSet(vec, 1.0);
  value = petsc_vec_to_adj_vector(&vec);
  ierr = adj_storage_memory_copy(value, &storage);
  for (i = 0; i < 5; i++)
  {
    ierr = adj_record_variable(&adjointer, u[i], storage);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");
  }

  /* One last spanner in the works. Claim that u20 is necessary for timestep 1's functional computation. */
  ierr = adj_timestep_set_functional_dependencies(&adjointer, 1, "Drag", 1, &u[3]);
  adj_test_assert(ierr == ADJ_OK, "Should have worked fine");

  has_val = has_values(&adjointer, 5, u);
  adj_test_assert(memcmp(has_val, before_2, 5 * sizeof(int)) == 0, "Should be {1, 1, 1, 1, 1}");
  free(has_val);

  ierr = adj_forget_adjoint_equation(&adjointer, 5);
  adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "Should not have worked");

  ierr = adj_forget_adjoint_equation(&adjointer, 4);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  has_val = has_values(&adjointer, 5, u);
  adj_test_assert(memcmp(has_val, before_2, 5 * sizeof(int)) == 0, "Should be {1, 1, 1, 1, 1}");
  free(has_val);

  ierr = adj_forget_adjoint_equation(&adjointer, 3);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  has_val = has_values(&adjointer, 5, u);
  adj_test_assert(memcmp(has_val, before_2, 5 * sizeof(int)) == 0, "Should be {1, 1, 1, 1, 1}");
  free(has_val);

  ierr = adj_forget_adjoint_equation(&adjointer, 2);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  has_val = has_values(&adjointer, 5, u);
  adj_test_assert(memcmp(has_val, after_2, 5 * sizeof(int)) == 0, "Should be {1, 1, 1, 1, 0}");
  free(has_val);

  ierr = adj_forget_adjoint_equation(&adjointer, 1);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  has_val = has_values(&adjointer, 5, u);
  adj_test_assert(memcmp(has_val, after_1, 5 * sizeof(int)) == 0, "Should be {1, 1, 1, 0, 0}");
  free(has_val);

  ierr = adj_forget_adjoint_equation(&adjointer, 0);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  has_val = has_values(&adjointer, 5, u);
  adj_test_assert(memcmp(has_val, after_0, 5 * sizeof(int)) == 0, "Should be {0, 0, 0, 0, 0}");
  free(has_val);
}

int* has_values(adj_adjointer *adjointer, int nb_vars, adj_variable *vars)
{
  int ierr, i = 0;
  int* ret = (int*) malloc(nb_vars * sizeof(int));
  adj_variable_data* data_ptr;
  for (i = 0;  i < nb_vars; i++)
  {
    ierr = adj_find_variable_data(&(adjointer->varhash), &vars[i], &data_ptr);
    if (ierr != ADJ_OK) adj_test_assert(ierr == ADJ_OK, "Should have passed");
    ret[i] = data_ptr->storage.has_value;
  }

  return ret;
}



#endif
