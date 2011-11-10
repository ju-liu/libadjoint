#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_debug.h"
#include "libadjoint/adj_test_main.h"
#include <string.h>

#ifndef HAVE_PETSC
void test_adj_test_nonlinear_derivative_action_transpose(void)
{
  adj_test_assert(1 == 1, "Don't have PETSc so can't run this test.");
}
#else
#include "libadjoint/adj_petsc_data_structures.h"
#include "libadjoint/adj_petsc.h"

void nonlinear_derivative_action_callback(int ndepends, adj_variable* variables, adj_vector* dependencies, adj_variable derivative, adj_vector contraction, int hermitian, adj_vector input, adj_scalar coefficient, void* context, adj_vector* output);

void test_adj_test_nonlinear_derivative_action_transpose(void)
{
  adj_adjointer adjointer;
  adj_nonlinear_block nblock;
  adj_nonlinear_block_derivative nblock_deriv;
  int ierr, dim=10, number_of_tests = 1;
  Vec contraction_vec, input, output, u0_vec;
  adj_scalar tol = 1e-10;
  adj_vector contraction, u0;
  adj_variable var1, var2;
  adj_storage_data storage;

  adj_create_adjointer(&adjointer);
  adj_set_petsc_data_callbacks(&adjointer);
  VecCreateSeq(PETSC_COMM_SELF, dim, &contraction_vec);
  VecCreateSeq(PETSC_COMM_SELF, dim, &input);
  VecCreateSeq(PETSC_COMM_SELF, dim, &output);
  contraction = petsc_vec_to_adj_vector(&contraction_vec);
  petsc_vec_set_random_proc(&contraction);
  ierr = adj_register_operator_callback(&adjointer, ADJ_NBLOCK_DERIVATIVE_ACTION_CB, "NonlinearOperator", (void (*)(void)) nonlinear_derivative_action_callback);
  adj_test_assert(ierr==ADJ_OK, "Should have worked");
  /* Create a initial velocity variable and record it */
  ierr = adj_create_variable("Velocity", 0, 0, 0, &var1);
  adj_test_assert(ierr==ADJ_OK, "Should have worked");
  VecCreateSeq(PETSC_COMM_SELF, dim, &u0_vec);
  u0 = petsc_vec_to_adj_vector(&u0_vec);
  petsc_vec_set_random_proc(&u0);
  ierr = adj_storage_memory_copy(u0, &storage);
  adj_test_assert(ierr==ADJ_OK, "Should have worked");
  ierr = adj_record_variable(&adjointer, var1, storage);
  adj_test_assert(ierr==ADJ_OK, "Should have worked");
  ierr = adj_create_variable("Velocity", 1, 0, 0, &var2);
  adj_test_assert(ierr==ADJ_OK, "Should have worked");
  ierr = adj_create_nonlinear_block("NonlinearOperator", 1, &var1, NULL, 1.0, &nblock);
  adj_test_assert(ierr==ADJ_OK, "Should have worked");
  ierr = adj_nonlinear_block_set_test_hermitian(&nblock, ADJ_TRUE, number_of_tests, tol); 
  adj_test_assert(ierr==ADJ_OK, "Should have worked");
  ierr = adj_create_nonlinear_block_derivative(&adjointer, nblock, var2, contraction, ADJ_FALSE, &nblock_deriv);
  adj_test_assert(ierr==ADJ_OK, "Should have worked");

  ierr = adj_test_nonlinear_derivative_action_transpose(&adjointer, nblock_deriv, petsc_vec_to_adj_vector(&input), petsc_vec_to_adj_vector(&output), number_of_tests, tol);
  adj_chkierr(ierr);
}

void nonlinear_derivative_action_callback(int ndepends, adj_variable* variables, adj_vector* dependencies, adj_variable derivative, adj_vector contraction, int hermitian, adj_vector input, adj_scalar coefficient, void* context, adj_vector* output)
{
  (void) hermitian;
  (void) context;
  (void) ndepends;
  (void) variables;
  (void) dependencies;
  (void) derivative;
  Vec input_shifted, input_vec;
  int i, shift = 0, dim;
  PetscScalar  *input_array, *shifted_array;
  Vec *output_vec;
  output_vec = (Vec*) malloc(sizeof(Vec));

  /* Apply the contraction in the non hermitian case */
  /* The contraction is simply a pointwise multiplication with the contraction vector */ 
  if (hermitian) 
  {
    VecDuplicate(petsc_vec_from_adj_vector(input), &input_vec);
    VecCopy(petsc_vec_from_adj_vector(input), input_vec);
  }
  else
  {
    VecDuplicate(petsc_vec_from_adj_vector(input), &input_vec);
    VecPointwiseMult(input_vec, petsc_vec_from_adj_vector(input), petsc_vec_from_adj_vector(contraction));
  }
  /* This is the identity operator */
  VecDuplicate(input_vec, output_vec);
  VecCopy(input_vec, *output_vec);
  /* Now, let us add an off diagonal term */
  VecDuplicate(input_vec, &input_shifted);
  VecCopy(input_vec, input_shifted);
  VecGetLocalSize(input_shifted, &dim);
  if (!hermitian) 
    shift = 2;
  else 
    shift = dim - 2;
  VecGetArray(input_shifted, &shifted_array);
  VecGetArray(input_vec, &input_array);
  for (i=0; i<dim; i++) {
      shifted_array[i] = input_array[(i+shift)%dim]; 
  }
  VecRestoreArray(input_shifted, &shifted_array);
 
  VecAXPY(*output_vec, 0.5, input_shifted);
  /* Apply the contraction in the hermitian case */
  if (hermitian) 
    VecPointwiseMult(*output_vec, *output_vec, petsc_vec_from_adj_vector(contraction));

  VecScale(*output_vec, (PetscScalar) coefficient);
  *output = petsc_vec_to_adj_vector(output_vec);
  
}
#endif
