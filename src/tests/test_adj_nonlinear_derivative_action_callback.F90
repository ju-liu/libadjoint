subroutine block_action_callback(ndepends, variables, dependencies, hermitian, coefficient, input, context, output_vec) bind(c)
  use libadjoint
  use libadjoint_petsc_data_structures
  use iso_c_binding
  implicit none
#include "libadjoint/adj_fortran.h"
#include "libadjoint/adj_petsc_f.h"

  integer(kind=c_int), intent(in), value :: ndepends
  type(adj_variable), dimension(ndepends), intent(in) :: variables
  type(adj_vector), dimension(ndepends), intent(in) :: dependencies
  integer(kind=c_int), intent(in), value :: hermitian
  adj_scalar_f, intent(in), value :: coefficient
  type(adj_vector), intent(in), value :: input
  type(c_ptr), intent(in), value :: context
  type(adj_vector), intent(out) :: output_vec

  integer, parameter :: m = 10
  integer :: ierr
  
#ifdef HAVE_PETSC
  Vec :: output
  call VecDuplicate(petsc_vec_from_adj_vector(input), output, ierr)
  call VecPointwiseMult(output, petsc_vec_from_adj_vector(input), petsc_vec_from_adj_vector(dependencies(1)), ierr)
  output_vec = petsc_vec_to_adj_vector(output)
#endif

end subroutine block_action_callback

subroutine nonlinear_derivative_action_callback(ndepends, variables, dependencies, derivative, contraction, hermitian, &
                                                & input, coefficient, context, output_vec) bind(c)
  use libadjoint
  use libadjoint_petsc_data_structures
  use iso_c_binding
  implicit none
#include "libadjoint/adj_fortran.h"
#include "libadjoint/adj_petsc_f.h"

  integer(kind=c_int), intent(in), value :: ndepends
  type(adj_variable), dimension(ndepends), intent(in) :: variables
  type(adj_vector), dimension(ndepends), intent(in) :: dependencies
  type(adj_variable), intent(in), value :: derivative
  type(adj_vector), intent(in), value :: contraction
  integer(kind=c_int), intent(in), value :: hermitian
  type(adj_vector), intent(in), value :: input
  adj_scalar_f, intent(in), value :: coefficient
  type(c_ptr), intent(in), value :: context
  type(adj_vector), intent(out) :: output_vec

  integer, parameter :: m = 10
  integer :: ierr

#ifdef HAVE_PETSC
  Vec :: output
  call VecDuplicate(petsc_vec_from_adj_vector(input), output, ierr)
  call VecPointwiseMult(output, petsc_vec_from_adj_vector(input), petsc_vec_from_adj_vector(dependencies(1)), ierr)
  output_vec = petsc_vec_to_adj_vector(output)
#endif
end subroutine nonlinear_derivative_action_callback

subroutine identity_assembly_callback(ndepends, variables, dependencies, hermitian, coefficient, & 
                                    & context, output_mat, rhs_vec) bind(c)
  use libadjoint
  use libadjoint_petsc_data_structures
  use iso_c_binding
  implicit none
#include "libadjoint/adj_fortran.h"
#include "libadjoint/adj_petsc_f.h"

  integer(kind=c_int), intent(in), value :: ndepends
  type(adj_variable), dimension(ndepends), intent(in) :: variables
  type(adj_vector), dimension(ndepends), intent(in) :: dependencies
  integer(kind=c_int), intent(in), value :: hermitian
  adj_scalar_f, intent(in), value :: coefficient
  type(c_ptr), intent(in), value :: context
  type(adj_matrix), intent(out) :: output_mat
  type(adj_vector), intent(out) :: rhs_vec

  integer, parameter :: m = 10
  integer :: ierr
  adj_scalar_f, parameter :: one = 1.0

#ifdef HAVE_PETSC
  Mat :: output
  Vec :: rhs, ones
  call MatCreateSeqAIJ(PETSC_COMM_SELF, m, m, 1, PETSC_NULL_INTEGER, output, ierr)
  call VecCreateSeq(PETSC_COMM_SELF, m, ones, ierr)
  call VecCreateSeq(PETSC_COMM_SELF, m, rhs, ierr)
  call VecZeroEntries(rhs, ierr)
  call VecSet(ones, one, ierr)
  call MatDiagonalSet(output, ones, INSERT_VALUES, ierr)
  call VecDestroy(ones, ierr)
  if (hermitian == ADJ_TRUE) then
    call MatHermitianTranspose(output, MAT_REUSE_MATRIX, output, ierr)
  end if
  rhs_vec = petsc_vec_to_adj_vector(rhs)
  output_mat = petsc_mat_to_adj_matrix(output)
#endif
end subroutine identity_assembly_callback

subroutine drag_derivative_callback(adjointer, variable, ndepends, dependencies, values, name, output) bind(c)
  use iso_c_binding
  use libadjoint_data_structures
  use libadjoint_petsc_data_structures
  implicit none
#include "libadjoint/adj_fortran.h"
#include "libadjoint/adj_petsc_f.h"

  type(adj_adjointer), intent(in) :: adjointer
  type(adj_variable), intent(in), value :: variable
  integer(kind=c_int), intent(in), value :: ndepends
  type(adj_variable), dimension(ndepends), intent(in) :: dependencies
  type(adj_vector), dimension(ndepends), intent(in) :: values
  character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
  type(adj_vector), intent(out) :: output
  integer, parameter :: m = 10
#ifdef HAVE_PETSC
  Vec :: output_petsc
  integer :: ierr

  call VecCreateSeq(PETSC_COMM_SELF, m, output_petsc, ierr)
  call VecZeroEntries(output_petsc, ierr)
  output = petsc_vec_to_adj_vector(output_petsc)
#endif
end subroutine drag_derivative_callback

subroutine test_adj_nonlinear_derivative_action_callback
  use libadjoint
  use libadjoint_petsc_data_structures
  use libadjoint_data_structures
  implicit none
#include "libadjoint/adj_fortran.h"
#include "libadjoint/adj_petsc_f.h"

  type(adj_adjointer) :: adjointer
  type(adj_block) :: N, I
  type(adj_nonlinear_block) :: N_nl
  type(adj_variable) :: u0, u1, adj_var0, adj_var1
  type(adj_equation) :: equation
  integer :: ierr
  integer, parameter :: m = 10
  type(adj_matrix) :: lhs
  type(adj_vector) :: rhs
  type(adj_vector) :: u0_vec, ones_vec
  type(adj_storage_data) :: storage
  procedure(adj_block_action_proc), bind(c) :: block_action_callback
  procedure(adj_nonlinear_derivative_action_proc), bind(c) :: nonlinear_derivative_action_callback
  procedure(adj_block_assembly_proc), bind(c) :: identity_assembly_callback
  procedure(adj_functional_derivative_proc), bind(c) :: drag_derivative_callback
  adj_scalar_f :: norm
  adj_scalar_f, parameter :: two = 2.0
  adj_scalar_f, parameter :: one = 1.0
#ifdef HAVE_PETSC
  Vec :: ones, u0_val
  PetscRandom :: rctx
#endif
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! In this example, the forward system is:
!! Au= (I         0) (u_0)
!!     (diag(u0)  I) (u_1)
!!
!! I.e. we have:
!! G= (0        0)
!!    (diag(u0) 0)
!!
!! And therefore the adjoint system reads as:
!! (A+G)*=(I 2diag(u0))
!!        (0 I        )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ierr = adj_create_adjointer(adjointer)
  call adj_chkierr(ierr)
  ierr = adj_set_petsc_data_callbacks(adjointer)

  ierr = adj_create_variable("Velocity", 0, 0, .false., u0)
  call adj_chkierr(ierr)
  ierr = adj_create_block("IdentityOperator", block=I)
  call adj_chkierr(ierr)
  ierr = adj_create_equation(u0, (/I/), (/u0/), equation)
  call adj_chkierr(ierr)
  ierr = adj_register_equation(adjointer, equation)
  call adj_chkierr(ierr)
  ierr = adj_destroy_equation(equation)
  call adj_chkierr(ierr)

  ierr = adj_create_nonlinear_block("NonlinearOperator", depends=(/u0/), nblock=N_nl)
  call adj_chkierr(ierr)
  ierr = adj_create_block("NonlinearOperator", nblock=N_nl, block=N)
  call adj_chkierr(ierr)
  ierr = adj_create_variable("Velocity", 1, 0, .false., u1)
  call adj_chkierr(ierr)
  ierr = adj_create_equation(u1, (/N, I/), (/u0, u1/), equation)
  call adj_chkierr(ierr)
  ierr = adj_register_equation(adjointer, equation)
  call adj_chkierr(ierr)
  ierr = adj_destroy_equation(equation)
  call adj_chkierr(ierr)

  ierr = adj_register_functional_derivative_callback(adjointer, "Drag", c_funloc(drag_derivative_callback))
  call adj_chkierr(ierr)
  ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, "IdentityOperator", &
  & c_funloc(identity_assembly_callback))
  ierr = adj_get_adjoint_equation(adjointer, functional="Drag", equation=1, lhs=lhs, rhs=rhs, adj_var=adj_var1)
  call adj_test_assert(ierr == ADJ_WARN_UNINITIALISED_VALUE, "Expected ADJ_WARN_UNINITIALISED_VALUE because no dependencies have been recorded.")

  ! Set adj_var1 to be constant 1.
#ifdef HAVE_PETSC
  call VecCreateSeq(PETSC_COMM_SELF, m, ones, ierr)
  call VecSet(ones, one, ierr)
#endif
  ones_vec = petsc_vec_to_adj_vector(ones)
  ierr = adj_storage_memory_copy(ones_vec, storage)
  ierr = adj_record_variable(adjointer, adj_var1, storage)
  call adj_test_assert(ierr == 0, "Should have worked")
 
  ! Lets create a random solution vector u1_val and record it
#ifdef HAVE_PETSC
  call VecCreateSeq(PETSC_COMM_SELF, m, u0_val, ierr)
  call PetscRandomCreate(PETSC_COMM_SELF, rctx, ierr)
  call PetscRandomSetFromOptions(rctx, ierr)
  call VecSetRandom(u0_val, rctx, ierr)
  call PetscRandomDestroy(rctx, ierr)
#endif
  u0_vec = petsc_vec_to_adj_vector(u0_val)
  ierr = adj_storage_memory_copy(u0_vec, storage)
  ierr = adj_record_variable(adjointer, u0, storage)
  call adj_test_assert(ierr == 0, "Should have worked")

  ierr = adj_get_adjoint_equation(adjointer, functional="Drag", equation=0, lhs=lhs, rhs=rhs, adj_var=adj_var0)
  call adj_test_assert(ierr == ADJ_ERR_NEED_CALLBACK, "Should not have worked")

  ! Register the callback for the block action
  ierr = adj_register_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, "NonlinearOperator", &
         & c_funloc(block_action_callback))
  call adj_test_assert(ierr == ADJ_OK, "Should have worked")

  ierr = adj_get_adjoint_equation(adjointer, functional="Drag", equation=0, lhs=lhs, rhs=rhs, adj_var=adj_var0)
  call adj_test_assert(ierr /= ADJ_OK, "Should not have worked")

  ! Register the callback for the nonlinear derivative action
  ierr = adj_register_operator_callback(adjointer, ADJ_NBLOCK_DERIVATIVE_ACTION_CB, "NonlinearOperator", &
                           & c_funloc(nonlinear_derivative_action_callback))
  call adj_test_assert(ierr == ADJ_OK, "Should have worked")

  ierr = adj_get_adjoint_equation(adjointer, functional="Drag", equation=0, lhs=lhs, rhs=rhs, adj_var=adj_var0)
  call adj_test_assert(ierr == ADJ_WARN_UNINITIALISED_VALUE, "Expected ADJ_WARN_UNINITIALISED_VALUE because no dependencies have been recorded.")
 
   ! The rhs should now be -2*u0_val
  call VecAXPY(petsc_vec_from_adj_vector(rhs), two, u0_val, ierr)
  call VecNorm(petsc_vec_from_adj_vector(rhs), NORM_2, norm, ierr)
  call adj_test_assert(norm==0.0, "The error should be zero!")

  ! Tidy up
#ifdef HAVE_PETSC
  call VecDestroy(ones, ierr)
  call VecDestroy(u0_val, ierr)
  call VecDestroy(petsc_vec_from_adj_vector(rhs), ierr)
  call MatDestroy(petsc_mat_from_adj_matrix(lhs), ierr)
#endif
  ierr = adj_destroy_adjointer(adjointer)

end subroutine test_adj_nonlinear_derivative_action_callback

