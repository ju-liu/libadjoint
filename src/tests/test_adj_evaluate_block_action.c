#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_evaluation.h"
#include "libadjoint/adj_test_main.h"
#include <string.h>

#ifndef HAVE_PETSC
void test_adj_evaluate_block_action(void)
{
  adj_test_assert(1 == 1, "Don't have PETSc so can't run this test.");
}
#else
#include "libadjoint/adj_petsc_data_structures.h"
#include "petsc.h"

void identity_action_callback(int nb_variables, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_vector input, void* context, adj_vector output);

void test_adj_evaluate_block_action(void)
{
  adj_adjointer adjointer;
  adj_block I;
  int ierr, dim=10;
  PetscReal norm;
  Vec input, output, difference;

  adj_create_adjointer(&adjointer);
  adj_set_petsc_data_callbacks(&adjointer);
  adj_register_operator_callback(&adjointer, ADJ_BLOCK_ACTION_CB, "IdentityOperator", (void (*)(void)) identity_action_callback);

  adj_create_block("IdentityOperator", NULL, NULL, &I);
  VecCreateSeq(PETSC_COMM_SELF, dim, &input);
  VecSet(input, 1.0);

  strncpy(I.name, "NotTheIdentityOperator", 22);
  I.name[22]='\0';
  ierr = adj_evaluate_block_action(&adjointer, I, petsc_vec_to_adj_vector(&input), (adj_vector*) &output);
  adj_test_assert(ierr!=ADJ_ERR_OK, "Should have not worked");

  strncpy(I.name, "IdentityOperator", 19);
  I.name[19]='\0';
  ierr = adj_evaluate_block_action(&adjointer, I, petsc_vec_to_adj_vector(&input), (adj_vector*) &output);
  adj_test_assert(ierr==ADJ_ERR_OK, "Should have worked");
 
  VecDuplicate(input, &difference);
  VecCopy(input, difference);
  VecAXPY(difference, -1.0, output);
  VecNorm(difference, NORM_2, &norm);
  adj_test_assert(norm == 0.0, "Norm should be zero");
}

void identity_action_callback(int nb_variables, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_vector input, void* context, adj_vector output)
{
  (void) hermitian;
  (void) context;
  (void) nb_variables;
  (void) variables;
  (void) dependencies;
  VecDuplicate(petsc_vec_from_adj_vector(input), (Vec*)output.ptr);
  VecCopy(petsc_vec_from_adj_vector(input), *(Vec*)output.ptr);
}
#endif
