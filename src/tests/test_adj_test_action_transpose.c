#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_debug.h"
#include "libadjoint/adj_test_main.h"
#include <string.h>

#ifndef HAVE_PETSC
void test_adj_evaluate_block_action(void)
{
  adj_test_assert(1 == 1, "Don't have PETSc so can't run this test.");
}
#else
#include "libadjoint/adj_petsc_data_structures.h"
#include "libadjoint/adj_petsc.h"


void matrix_action_callback(int nb_variables, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_scalar coefficient, adj_vector input, void* context, adj_vector* output);

void test_adj_test_action_transpose(void)
{
  adj_adjointer adjointer;
  adj_block matrix;
  int ierr, dim=10, number_of_tests = 1;
  Vec input, output;
  adj_scalar tol = 1e-10;

  adj_create_adjointer(&adjointer);
  adj_set_petsc_data_callbacks(&adjointer);
  adj_register_operator_callback(&adjointer, ADJ_BLOCK_ACTION_CB, "MatrixOperator", (void (*)(void)) matrix_action_callback);

  adj_create_block("MatrixOperator", NULL, NULL, &matrix);
  VecCreateSeq(PETSC_COMM_SELF, dim, &input);
  VecCreateSeq(PETSC_COMM_SELF, dim, &output);
  VecSet(input, 1.0);
  VecSet(output, 1.0);
 
  ierr = adj_test_action_transpose(&adjointer, matrix, petsc_vec_to_adj_vector(&input), petsc_vec_to_adj_vector(&output), number_of_tests, tol);
  adj_test_assert(ierr==ADJ_ERR_OK, "Should have worked");
}

void matrix_action_callback(int nb_variables, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_scalar coefficient, adj_vector input, void* context, adj_vector* output)
{
  (void) hermitian;
  (void) context;
  (void) nb_variables;
  (void) variables;
  (void) dependencies;
  VecDuplicate(petsc_vec_from_adj_vector(input), (Vec*)output->ptr);
  VecCopy(petsc_vec_from_adj_vector(input), *(Vec*)output->ptr);
  VecScale(*(Vec*)output->ptr, (PetscScalar) coefficient);
}
#endif
