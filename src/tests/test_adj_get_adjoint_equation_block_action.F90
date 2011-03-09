#ifndef HAVE_PETSC
subroutine test_adj_get_adjoint_equation_block_action

  use libadjoint
  use iso_c_binding
  implicit none

  call adj_test_assert(1 == 1, "Don't have PETSc so can't run this test")
end subroutine test_adj_get_adjoint_equation_block_action
#else

subroutine identity_assembly_callback(nvar, variables, dependencies, hermitian, coefficient, context, output, rhs) bind(c)
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
  adj_scalar_f, intent(in), value :: coefficient
  type(c_ptr), intent(in), value :: context
  type(adj_matrix), intent(out) :: output
  type(adj_vector), intent(out) :: rhs

  integer, parameter :: m = 2
  integer :: ierr

  PetscScalar, parameter :: one = 1.0
  PetscScalar :: coefficient_petsc
  Vec :: ones
  Mat :: output_mat
  Vec :: rhs_petsc

  call adj_test_assert(nvar == 0, "We don't depend on any variables")

  call MatCreateSeqAIJ(PETSC_COMM_SELF, m, m, 1, PETSC_NULL_INTEGER, output_mat, ierr)
  call MatZeroEntries(output_mat, ierr)
  call VecCreateSeq(PETSC_COMM_SELF, m, ones, ierr)
  call VecSet(ones, one, ierr)
  call MatDiagonalSet(output_mat, ones, INSERT_VALUES, ierr)
  call VecDestroy(ones, ierr)
! MatHermitianTranspose was only added in petsc 3.1, it seems
#if PETSC_VERSION_MINOR > 0
  if (hermitian == ADJ_TRUE) then
    call MatHermitianTranspose(output_mat, MAT_REUSE_MATRIX, output_mat, ierr)
  end if
#endif
  coefficient_petsc = coefficient ! I wish fortran had casts!
  call MatScale(output_mat, coefficient_petsc, ierr)
  output = petsc_mat_to_adj_matrix(output_mat)
  call MatGetVecs(output_mat, PETSC_NULL_OBJECT, rhs_petsc, ierr)
  call VecZeroEntries(rhs_petsc, ierr)
  rhs = petsc_vec_to_adj_vector(rhs_petsc)
end subroutine identity_assembly_callback

subroutine identity_action_callback(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
  use iso_c_binding
  use libadjoint
  use libadjoint_petsc_data_structures
  implicit none
#include "libadjoint/adj_petsc_f.h"

  integer(kind=c_int), intent(in), value :: nvar
  type(adj_variable), dimension(nvar), intent(in) :: variables
  type(adj_vector), dimension(nvar), intent(in) :: dependencies
  integer(kind=c_int), intent(in), value :: hermitian
  adj_scalar_f, intent(in), value :: coefficient
  type(adj_vector), intent(in), value :: input
  type(c_ptr), intent(in), value :: context
  type(adj_vector), intent(out) :: output

  Vec :: output_vec
  integer :: ierr
  PetscScalar :: coefficient_petsc

  call adj_test_assert(nvar == 0, "We don't depend on any variables")
  call VecDuplicate(petsc_vec_from_adj_vector(input), output_vec, ierr)
  call VecCopy(petsc_vec_from_adj_vector(input), output_vec, ierr)
  coefficient_petsc = coefficient ! I wish fortran had casts!
  call VecScale(output_vec, coefficient_petsc, ierr)
  output = petsc_vec_to_adj_vector(output_vec)
end subroutine identity_action_callback

subroutine functional_derivative_callback(variable, ndepends, dependencies, values, &
                                        & name, start_time, end_time, output) bind(c)
  use iso_c_binding
  use libadjoint
  use libadjoint_petsc_data_structures
  implicit none
#include "libadjoint/adj_petsc_f.h"

  type(adj_variable), intent(in), value :: variable
  integer(kind=c_int), intent(in), value :: ndepends
  type(adj_variable), dimension(ndepends), intent(in) :: dependencies
  type(adj_vector), dimension(ndepends), intent(in) :: values
  character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
  adj_scalar_f, intent(in), value :: start_time
  adj_scalar_f, intent(in), value :: end_time
  type(adj_vector), intent(out) :: output

  Vec :: output_vec
  integer, parameter :: m = 2
  integer :: ierr

  call adj_test_assert(ndepends == 1, "We depend on one variable")
  call VecCreateSeq(PETSC_COMM_SELF, m, output_vec, ierr)
  call VecZeroEntries(output_vec, ierr)
  output = petsc_vec_to_adj_vector(output_vec)

end subroutine functional_derivative_callback

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
  procedure(adj_functional_derivative_proc), bind(c) :: functional_derivative_callback

  PetscScalar, parameter :: one = 1.0
  type(adj_matrix) :: lhs
  type(adj_vector) :: rhs
  Vec, target :: lambda0, lambda1, sum, u0_vec
  type(adj_vector) :: lambda1_vec
  PetscScalar :: norm
  PetscRandom :: rctx
  KSP :: ksp

  ierr = adj_create_adjointer(adjointer)
  
  ! Test the html output
  ierr = adj_adjointer_to_html(adjointer, "test_adj_get_adjoint_equation_block_action_0_forward.html", ADJ_FORWARD)
  call adj_test_assert(ierr == ADJ_ERR_OK, "html output should have worked")
  ierr = adj_adjointer_to_html(adjointer, "test_adj_get_adjoint_equation_block_action_0_adjoint.html", ADJ_ADJOINT)
  call adj_test_assert(ierr == ADJ_ERR_OK, "html output should have worked")
  
  call adj_chkierr(ierr)
  ierr = adj_set_petsc_data_callbacks(adjointer)
  call adj_chkierr(ierr)
  ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "IdentityOperator", c_funloc(identity_assembly_callback))
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ierr = adj_create_block("IdentityOperator", block=I)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ierr = adj_create_variable("Velocity", timestep=0, iteration=0, auxiliary=ADJ_FALSE, variable=u0)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ierr = adj_create_variable("Velocity", timestep=1, iteration=0, auxiliary=ADJ_FALSE, variable=u1)
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

  ierr = adj_get_adjoint_equation(adjointer, equation=2, functional="Drag", lhs=lhs, rhs=rhs, variable=adj_var1)
  call adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "Should not have worked")

  ierr = adj_get_adjoint_equation(adjointer, equation=1, functional="Drag", lhs=lhs, rhs=rhs, variable=adj_var1)
  call adj_test_assert(ierr == ADJ_ERR_NEED_CALLBACK, "Need the functional callback")

  ierr = adj_register_functional_derivative_callback(adjointer, "Drag", c_funloc(functional_derivative_callback))
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ierr = adj_timestep_set_functional_dependencies(adjointer, timestep=0, functional="Drag", dependencies=(/u0/))
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")
  ierr = adj_timestep_set_functional_dependencies(adjointer, timestep=1, functional="Drag", dependencies=(/u0/))
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ierr = adj_get_adjoint_equation(adjointer, equation=1, functional="Drag", lhs=lhs, rhs=rhs, variable=adj_var1)
  call adj_test_assert(ierr == ADJ_ERR_NEED_VALUE, "We should need the value for u0")

  ierr = adj_timestep_set_functional_dependencies(adjointer, timestep=1, functional="Drag", dependencies=(/u0, u1/))
  call adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "We can't set the functional dependencies twice")

  ! Test the html output
  ierr = adj_adjointer_to_html(adjointer, "test_adj_get_adjoint_equation_block_action_1_forward.html", ADJ_FORWARD)
  call adj_test_assert(ierr == ADJ_ERR_OK, "html output should have worked")
  ierr = adj_adjointer_to_html(adjointer, "test_adj_get_adjoint_equation_block_action_1_adjoint.html", ADJ_ADJOINT)

  call adj_test_assert(ierr == ADJ_ERR_OK, "html output should have worked")

  ! Record a value for u0
  call VecCreateSeq(PETSC_COMM_SELF, m, u0_vec, ierr)
  call VecZeroEntries(u0_vec, ierr)
  ierr = adj_record_variable(adjointer, u0, adj_storage_memory_copy(petsc_vec_to_adj_vector(u0_vec)))
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ierr = adj_get_adjoint_equation(adjointer, equation=1, functional="Drag", lhs=lhs, rhs=rhs, variable=adj_var1)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ! We don't actually need the memory for lhs and rhs, so we'll delete them now
  call petsc_vec_destroy_proc(rhs)
  call petsc_mat_destroy_proc(lhs)

  ierr = adj_get_adjoint_equation(adjointer, equation=0, functional="Drag", lhs=lhs, rhs=rhs, variable=adj_var1)
  call adj_test_assert(ierr == ADJ_ERR_NEED_VALUE, "Should need the value for lambda1")

  ! We'll decide on a random value for lambda1 (lambda1 = dJ/du, so that's the same
  ! as a random functional)
  call VecCreateSeq(PETSC_COMM_SELF, m, lambda1, ierr)
  call PetscRandomCreate(PETSC_COMM_SELF, rctx, ierr)
  call PetscRandomSetFromOptions(rctx, ierr)
  call VecSetRandom(lambda1, rctx, ierr)
  call PetscRandomDestroy(rctx, ierr)

  lambda1_vec%klass = 0
  lambda1_vec%ptr = c_loc(lambda1)
  ierr = adj_record_variable(adjointer, adj_var1, adj_storage_memory_copy(lambda1_vec))
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ierr = adj_get_adjoint_equation(adjointer, equation=0, functional="Drag", lhs=lhs, rhs=rhs, variable=adj_var1)
  call adj_test_assert(ierr == ADJ_ERR_NEED_CALLBACK, "Should need the block_action_callback")

  ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "IdentityOperator", c_funloc(identity_action_callback))
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ierr = adj_get_adjoint_equation(adjointer, equation=0, functional="Drag", lhs=lhs, rhs=rhs, variable=adj_var1)
  call adj_test_assert(ierr == ADJ_ERR_OK, "Should have worked")

  ! So now solve lhs . lambda0 = rhs
  ! With this setup, lambda0 = -1 * lambda1.

  call VecDuplicate(lambda1, lambda0, ierr)
  call KSPCreate(PETSC_COMM_SELF, ksp, ierr)
  call KSPSetOperators(ksp, petsc_mat_from_adj_matrix(lhs), petsc_mat_from_adj_matrix(lhs), SAME_NONZERO_PATTERN, ierr)
  call KSPSetFromOptions(ksp, ierr)
  call KSPSolve(ksp, petsc_vec_from_adj_vector(rhs), lambda0, ierr)
  call petsc_vec_destroy_proc(rhs)
  call petsc_mat_destroy_proc(lhs)

  ! OK. Now compare lambda0 and lambda1.
  call VecDuplicate(lambda1, sum, ierr)
  call VecZeroEntries(sum, ierr)
  call VecAXPY(sum, one, lambda0, ierr)
  call VecAXPY(sum, one, lambda1, ierr)
  call VecNorm(sum, NORM_2, norm, ierr)

  call VecDestroy(lambda0, ierr)
  call VecDestroy(lambda1, ierr)
  call VecDestroy(sum, ierr)

  call adj_test_assert(norm == 0, "lambda1 should equal -1 * lambda0")
end subroutine test_adj_get_adjoint_equation_block_action
#endif
