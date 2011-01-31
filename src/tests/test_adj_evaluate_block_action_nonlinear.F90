#ifndef HAVE_PETSC
subroutine test_adj_evaluate_block_action_nonlinear

  use libadjoint
  use iso_c_binding
  implicit none

  call adj_test_assert(1 == 1, "Don't have PETSc so can't run this test.")
end subroutine test_adj_evaluate_block_action_nonlinear

#else

#include "libadjoint/adj_fortran.h"

subroutine timestepping_action_callback(nvar, variables, dependencies, hermitian, input, context, output) &
         & bind(c, name='timestepping_action_callback_')
  use libadjoint
  use iso_c_binding
  implicit none
#include "finclude/petsc.h"

  integer(kind=c_int), intent(in), value :: nvar
  type(adj_variable), dimension(nvar), intent(in) :: variables
  type(adj_vector), dimension(nvar), intent(in) :: dependencies
  integer(kind=c_int), intent(in), value :: hermitian
  type(adj_vector), intent(in), value :: input
  type(c_ptr), intent(in), value :: context
  type(adj_vector), intent(out) :: output

  integer :: ierr
  Vec, pointer :: input_vec, output_vec, dependency_vec

  call c_f_pointer(input%ptr, input_vec)
  call c_f_pointer(output%ptr, output_vec)
  call c_f_pointer(dependencies(1)%ptr, dependency_vec)

  call adj_test_assert(nvar == 1, "We only depend on one variable")

  call VecDuplicate(input_vec, output_vec, ierr)
  call VecPointwiseMult(output_vec, dependency_vec, input_vec, ierr)
end subroutine timestepping_action_callback

subroutine test_adj_evaluate_block_action_nonlinear
  use libadjoint
  use libadjoint_petsc_data_structures
  use iso_c_binding
  implicit none
#include "finclude/petsc.h"

  type(adj_adjointer) :: adjointer
  procedure(adj_block_action_proc) :: timestepping_action_callback
  integer :: ierr
  type(adj_block) :: T, I
  type(adj_nonlinear_block) :: V
  type(adj_variable) :: u0
  Vec :: input, output, difference
  Vec, target :: u0_val
  integer, parameter :: m = 2
  PetscScalar :: norm
  adj_scalar_f, dimension(m) :: values
  integer, dimension(m) :: indices
  integer :: j
  type(adj_equation) :: equation
  PetscScalar, parameter :: one = 1.0
  type(adj_vector) :: u0_vec

  ierr = adj_create_adjointer(adjointer)
  call adj_chkierr(ierr)
  ierr = adj_set_petsc_data_callbacks(adjointer)
  call adj_chkierr(ierr)
  ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, &
                                      & "TimesteppingOperator", c_funloc(timestepping_action_callback))
  call adj_chkierr(ierr)

  ierr = adj_create_variable("Velocity", 0, 0, ADJ_FALSE, u0)
  call adj_chkierr(ierr)
  ierr = adj_create_nonlinear_block("AdvectionOperator", depends=(/u0/), nblock=V)
  call adj_chkierr(ierr)
  ierr = adj_create_block("TimesteppingOperator", nblock=V, block=T)
  call adj_chkierr(ierr)
  ierr = adj_create_block("Identity", block=I)
  call adj_chkierr(ierr)
  ierr = adj_create_equation(var=u0, blocks=(/I/), targets=(/u0/), equation=equation)
  call adj_test_assert(ierr == 0, "Should have worked")

  ! Set up the input
  call VecCreateSeq(PETSC_COMM_SELF, m, input, ierr)
  call VecSet(input, one, ierr)

  !T%name = "NotTheTimesteppingOperator"
  !call adj_evaluate_block_action(adjointer, T, input, output, ierr)
  !call adj_test_assert(ierr == ADJ_ERR_NEED_CALLBACK, "Should need the callback")
  !T%name = "TimesteppingOperator"
  !call adj_evaluate_block_action(adjointer, T, input, output, ierr)
  !call adj_test_assert(ierr == ADJ_ERR_NEED_VALUE, "Should need the value of u0")

  ! Set up u0_val
  call VecCreateSeq(PETSC_COMM_SELF, m, u0_val, ierr)
  do j=1,m
    indices(j) = j-1
    values(j) = float(m)
  end do
  call VecSetValues(u0_val, m, indices, values, INSERT_VALUES, ierr)
  call VecAssemblyBegin(u0_val, ierr)
  call VecAssemblyEnd(u0_val, ierr)

  u0_vec%ptr = c_loc(u0_val)
  ierr = adj_record_variable(adjointer, u0, adj_storage_memory(u0_vec))
  call adj_test_assert(ierr == 0, "Should have worked")
  !call adj_evaluate_block_action(adjointer, T, input, output, ierr)
  call adj_test_assert(ierr == 0, "Should have worked")

  norm = 0.0
  call VecDuplicate(u0_val, difference, ierr)
  call VecCopy(u0_val, difference, ierr)
  call VecAXPY(difference, -one, output, ierr)
  call VecNorm(difference, NORM_2, norm, ierr)

  call adj_test_assert(norm == 0.0, "should be zero")

  ierr = adj_destroy_adjointer(adjointer)
  call VecDestroy(input, ierr)
  call VecDestroy(difference, ierr)
  call VecDestroy(output, ierr)
end subroutine test_adj_evaluate_block_action_nonlinear

#endif
