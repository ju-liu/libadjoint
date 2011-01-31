#ifndef HAVE_PETSC
subroutine test_adj_get_adjoint_equation_block_action

  use libadjoint
  use iso_c_binding
  implicit none

  call adj_test_assert(1 == 1, "Don't have PETSc so can't run this test")
end subroutine test_adj_get_adjoint_equation_block_action
#else

subroutine identity_assembly_callback(nvar, variables, dependencies, hermitian, context, output) bind(c)
  use iso_c_binding
  use libadjoint
  use libadjoint_petsc_data_structures
  implicit none
#include "libadjoint/adj_fortran.h"
#include "libadjoint/adj_petsc_f.h"

  integer(kind=c_int), intent(in), value :: nvar
  type(adj_variable), dimension(nvar), intent(in) :: variables
  type(adj_vector), dimension(nvar), intent(in) :: dependencies
  integer(kind=c_int), intent(in), value :: hermitian
  type(c_ptr), intent(in), value :: context
  type(adj_matrix), intent(out) :: output

  integer, parameter :: m = 2
  integer :: ierr

  PetscScalar, parameter :: one = 1.0
  Vec :: ones
  Mat :: output_mat

  call adj_test_assert(nvar == 0, "We don't depend on any variables")

  call MatCreateSeqAIJ(PETSC_COMM_SELF, m, m, 1, PETSC_NULL_INTEGER, output_mat, ierr)
  call MatZeroEntries(output_mat, ierr)
  call VecCreateSeq(PETSC_COMM_SELF, m, ones, ierr)
  call VecSet(ones, one, ierr)
  call MatDiagonalSet(output_mat, ones, INSERT_VALUES, ierr)
  call VecDestroy(ones, ierr)
  if (hermitian == ADJ_TRUE) then
    call MatHermitianTranspose(output_mat, MAT_REUSE_MATRIX, output_mat, ierr)
  end if
  output = petsc_mat_to_adj_matrix(output_mat)
end subroutine identity_assembly_callback

subroutine identity_action_callback(nvar, variables, dependencies, hermitian, input, context, output) bind(c)
  use iso_c_binding
  use libadjoint
  use libadjoint_petsc_data_structures
  implicit none
#include "libadjoint/adj_petsc_f.h"

  integer(kind=c_int), intent(in), value :: nvar
  type(adj_variable), dimension(nvar), intent(in) :: variables
  type(adj_vector), dimension(nvar), intent(in) :: dependencies
  integer(kind=c_int), intent(in), value :: hermitian
  type(adj_vector), intent(in), value :: input
  type(c_ptr), intent(in), value :: context
  type(adj_vector), intent(out) :: output

  Vec :: output_vec
  integer :: ierr

  call adj_test_assert(nvar == 0, "We don't depend on any variables")
  call VecDuplicate(petsc_vec_from_adj_vector(input), output_vec, ierr)
  call VecCopy(petsc_vec_from_adj_vector(input), output_vec, ierr)
  output = petsc_vec_to_adj_vector(output_vec)
end subroutine identity_action_callback

subroutine test_adj_get_adjoint_equation_block_action
  use libadjoint
  use libadjoint_petsc_data_structures
  implicit none
#include "libadjoint/adj_petsc_f.h"

  type(adj_adjointer) :: adjointer
  type(adj_block) :: I
  type(adj_variable) :: u0, u1, adj_var0, adj_var1
  type(adj_equation) :: equation
  integer :: ierr
  integer, parameter :: m = 2
  procedure(adj_block_assembly_proc), bind(c) :: identity_assembly_callback
  procedure(adj_block_action_proc), bind(c) :: identity_action_callback

  PetscScalar, parameter :: one = 1.0
  type(adj_matrix) :: lhs
  type(adj_vector) :: rhs
  Vec :: lambda0, lambda1, sum
  PetscScalar :: norm
  PetscRandom :: rctx
  KSP :: ksp

  ierr = adj_create_adjointer(adjointer)
  call adj_chkierr(ierr)
  ierr = adj_set_petsc_data_callbacks(adjointer)
  call adj_chkierr(ierr)
  ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "IdentityOperator", c_funloc(identity_assembly_callback))
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ierr = adj_create_block("IdentityOperator", block=I)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ierr = adj_create_variable("Velocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, var=u0)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ierr = adj_create_variable("Velocity", timestep=1, iteration=0, auxiliary=ADJ_FALSE, var=u1)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ierr = adj_create_equation(u0, (/I/), (/u0/), equation)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ierr = adj_register_equation(adjointer, equation)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ierr = adj_destroy_equation(equation)
  call adj_chkierr(ierr)

  ierr = adj_create_equation(u1, (/I, I/), (/u0, u1/), equation)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ierr = adj_register_equation(adjointer, equation)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ierr = adj_destroy_equation(equation)
  call adj_chkierr(ierr)

  ierr = adj_get_adjoint_equation(adjointer, equation=2, functional=0, lhs=lhs, rhs=rhs, var=adj_var1)
  call adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "Should not have worked")

  ierr = adj_get_adjoint_equation(adjointer, equation=1, functional=0, lhs=lhs, rhs=rhs, var=adj_var1)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ! We don't actually need the memory for lhs and rhs, so we'll delete them now
  call petsc_vec_destroy_proc(rhs)
  call petsc_mat_destroy_proc(lhs)

end subroutine test_adj_get_adjoint_equation_block_action
#endif
