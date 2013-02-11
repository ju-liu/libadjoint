#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_test_main.h"
#include "libadjoint/adj_adjointer_visualisation.h"

#ifndef HAVE_PETSC
void test_adj_variable_comparison(void)
{
  adj_test_assert(1 == 1, "Don't have PETSc so can't run this test.");
}
#else
#include "libadjoint/adj_petsc_data_structures.h"
void test_adj_variable_comparison(void)
{
  adj_adjointer adjointer;
  adj_variable u_var;
  adj_storage_data storage;
  int ierr;

  Vec u_petsc;

  adj_create_adjointer(&adjointer);
  adj_set_petsc_data_callbacks(&adjointer);

  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &u_var);

  VecCreateSeq(PETSC_COMM_SELF, 2, &u_petsc);
  VecSet(u_petsc, (PetscScalar) 1.0);

  ierr = adj_storage_memory_copy(petsc_vec_to_adj_vector(&u_petsc), &storage);
  adj_test_assert(ierr == ADJ_OK, "Should be OK");

  ierr = adj_record_variable(&adjointer, u_var, storage);
  adj_test_assert(ierr == ADJ_OK, "Should be OK");

  VecSet(u_petsc, (PetscScalar) 2.0);
  ierr = adj_storage_memory_copy(petsc_vec_to_adj_vector(&u_petsc), &storage);
  adj_test_assert(ierr == ADJ_OK, "Should be OK");

  ierr = adj_storage_set_compare(&storage, 999, 0.0);
  adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "Should warn about the 999");

  ierr = adj_storage_set_compare(&storage, ADJ_TRUE, -1.0);
  adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "Should warn about the -1.0");

  ierr = adj_storage_set_compare(&storage, ADJ_TRUE, 0.0);
  adj_test_assert(ierr == ADJ_OK, "Should be OK");

  ierr = adj_record_variable(&adjointer, u_var, storage);
  adj_test_assert(ierr == ADJ_WARN_COMPARISON_FAILED, "Should have failed the comparison");

  VecDestroy(&u_petsc);
  ierr = adj_destroy_adjointer(&adjointer);
}
#endif
