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
  int ierr, dim=10, number_of_tests = 5;
  Vec input, output;
  adj_scalar tol = 1e-10;

  adj_create_adjointer(&adjointer);
  adj_set_petsc_data_callbacks(&adjointer);
  adj_register_operator_callback(&adjointer, ADJ_BLOCK_ACTION_CB, "MatrixOperator", (void (*)(void)) matrix_action_callback);

  adj_create_block("MatrixOperator", NULL, NULL, &matrix);
  VecCreateSeq(PETSC_COMM_SELF, dim, &input);
  VecCreateSeq(PETSC_COMM_SELF, dim, &output);
 
  ierr = adj_test_block_action_transpose(&adjointer, matrix, petsc_vec_to_adj_vector(&input), petsc_vec_to_adj_vector(&output), number_of_tests, tol);
  adj_test_assert(ierr==ADJ_ERR_OK, "Should have worked");
}

void matrix_action_callback(int nb_variables, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_scalar coefficient, adj_vector input, void* context, adj_vector* output)
{
  (void) hermitian;
  (void) context;
  (void) nb_variables;
  (void) variables;
  (void) dependencies;
  Vec input_shifted;
  int i, shift = 0, dim;
  PetscScalar  *input_array, *shifted_array;
  Vec *output_vec;
  output_vec = (Vec*) malloc(sizeof(Vec));

  /* This is the identity operator */
  VecDuplicate(petsc_vec_from_adj_vector(input), output_vec);
  VecCopy(petsc_vec_from_adj_vector(input), *output_vec);
  /* Now, let us add an off diagonal term */
  VecDuplicate(petsc_vec_from_adj_vector(input), &input_shifted);
  VecCopy(petsc_vec_from_adj_vector(input), input_shifted);
  VecGetLocalSize(input_shifted, &dim);
  if (!hermitian) 
    shift = 2;
  else 
    shift = dim - 2;
  VecGetArray(input_shifted, &shifted_array);
  VecGetArray(petsc_vec_from_adj_vector(input), &input_array);
  for (i=0; i<dim; i++) {
      shifted_array[i] = input_array[(i+shift)%dim]; 
  }
  VecRestoreArray(input_shifted, &shifted_array);
 
  VecAXPY(*output_vec, 0.5, input_shifted);

  VecScale(*output_vec, (PetscScalar) coefficient);
  *output = petsc_vec_to_adj_vector(output_vec);
}
#endif
