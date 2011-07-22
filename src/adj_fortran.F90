#include "libadjoint/adj_fortran.h"
! Yes, writing this file really was as boring as you would imagine

module libadjoint_data_structures
  use iso_c_binding
  implicit none

  type, bind(c) :: adj_variable
    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name
    integer(kind=c_int) :: timestep
    integer(kind=c_int) :: iteration
    integer(kind=c_int) :: type
    integer(kind=c_int) :: auxiliary
    character(kind=c_char), dimension(ADJ_NAME_LEN) :: functional
  end type adj_variable

  type, bind(c) :: adj_nonlinear_block
    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name
    adj_scalar_f :: coefficient
    integer(kind=c_int) :: ndepends
    type(c_ptr) :: depends
    type(c_ptr) :: context
    integer(kind=c_int) :: test_deriv_hermitian
    integer(kind=c_int) :: number_of_tests
    adj_scalar_f :: tolerance
    integer(kind=c_int) :: test_derivative
    integer(kind=c_int) :: number_of_rounds
  end type adj_nonlinear_block

  type, bind(c) :: adj_block
    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name
    integer(kind=c_int) :: has_nonlinear_block
    type(adj_nonlinear_block) :: nonlinear_block
    type(c_ptr) :: context
    integer(kind=c_int) :: hermitian
    adj_scalar_f :: coefficient
    integer(kind=c_int) :: test_hermitian
    integer(kind=c_int) :: number_of_tests
    adj_scalar_f :: tolerance
  end type adj_block

  type, bind(c) :: adj_equation
    type(adj_variable) :: variable
    integer(kind=c_int) :: nblocks
    type(c_ptr) :: blocks
    type(c_ptr) :: targets
    integer(kind=c_int) :: nrhsdeps
    type(c_ptr) :: rhsdeps
    type(c_ptr) :: rhs_context
  end type adj_equation

  type, bind(c) :: adj_data_callbacks
    type(c_funptr) :: vec_duplicate
    type(c_funptr) :: vec_axpy
    type(c_funptr) :: vec_destroy
    type(c_funptr) :: vec_set_values
    type(c_funptr) :: vec_get_size
    type(c_funptr) :: vec_divide
    type(c_funptr) :: vec_get_norm
    type(c_funptr) :: vec_dot_product
    type(c_funptr) :: vec_set_random

    type(c_funptr) :: mat_duplicate
    type(c_funptr) :: mat_axpy
    type(c_funptr) :: mat_destroy
  end type adj_data_callbacks

  type, bind(c) :: adj_variable_data_list
    type(c_ptr) :: firstnode
    type(c_ptr) :: lastnode
  end type adj_variable_data_list

  type, bind(c) :: adj_op_callback_list
    type(c_ptr) :: firstnode
    type(c_ptr) :: lastnode
  end type adj_op_callback_list
  
  type, bind(c) :: adj_func_callback_list
    type(c_ptr) :: firstnode
    type(c_ptr) :: lastnode
  end type adj_func_callback_list

  type, bind(c) :: adj_func_deriv_callback_list
    type(c_ptr) :: firstnode
    type(c_ptr) :: lastnode
  end type adj_func_deriv_callback_list

  type, bind(c) :: adj_adjointer
    integer(kind=c_int) :: nequations
    integer(kind=c_int) :: equations_sz
    type(c_ptr) :: equations

    integer(kind=c_int) :: ntimesteps
    type(c_ptr) :: timestep_data

    type(c_ptr) :: varhash
    type(adj_variable_data_list) :: vardata

    integer(kind=c_int), dimension(ADJ_NO_OPTIONS) :: options

    type(adj_data_callbacks) :: callbacks
    type(adj_op_callback_list) :: nonlinear_colouring_list
    type(adj_op_callback_list) :: nonlinear_action_list
    type(adj_op_callback_list) :: nonlinear_derivative_action_list
    type(adj_op_callback_list) :: nonlinear_derivative_assembly_list
    type(adj_op_callback_list) :: block_action_list
    type(adj_op_callback_list) :: block_assembly_list
    type(adj_func_callback_list) :: functional_list
    type(adj_func_deriv_callback_list) :: functional_derivative_list
    type(c_funptr) :: forward_source_callback
  end type adj_adjointer

  type, bind(c) :: adj_vector
    type(c_ptr) :: ptr
    integer(kind=c_int) :: klass
    integer(kind=c_int) :: flags
  end type adj_vector

  type, bind(c) :: adj_matrix
    type(c_ptr) :: ptr
    integer(kind=c_int) :: klass
    integer(kind=c_int) :: flags
  end type adj_matrix

  type, bind(c) :: adj_storage_data
    integer(kind=c_int) :: storage_type
    integer(kind=c_int) :: has_value
    integer(kind=c_int) :: compare
    adj_scalar_f :: comparison_tolerance
    integer(kind=c_int) :: overwrite
    type(adj_vector) :: value
    type(c_ptr) :: filename
  end type adj_storage_data

  type, bind(c) :: adj_dictionary
    type(c_ptr) :: dict = c_null_ptr
  end type adj_dictionary
end module libadjoint_data_structures

module libadjoint

  use libadjoint_data_structures
  use iso_c_binding
  implicit none

  abstract interface
    subroutine adj_vec_duplicate_proc(x, newx) bind(c)
      ! Creates a new vector of the same type as an existing vector and set its entries to zero.
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_vector), intent(in), value :: x
      type(adj_vector), intent(out) :: newx
    end subroutine adj_vec_duplicate_proc

    subroutine adj_vec_axpy_proc(y, alpha, x) bind(c)
      ! Computes y = alpha x + y.
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_vector), intent(inout) :: y
      adj_scalar_f, intent(in), value :: alpha
      type(adj_vector), intent(in), value :: x
    end subroutine adj_vec_axpy_proc

    subroutine adj_vec_destroy_proc(x) bind(c)
      ! Destroys a vector.
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_vector), intent(inout) :: x
    end subroutine adj_vec_destroy_proc

    subroutine adj_vec_set_values_proc(vec, scalars) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_vector), intent(inout) :: vec
      adj_scalar_f, dimension(*), intent(in) :: scalars
    end subroutine adj_vec_set_values_proc

    subroutine adj_vec_pointwisedivide_proc(numerator, denominator) bind(c)
      use libadjoint_data_structures
      type(adj_vector), intent(inout) :: numerator
      type(adj_vector), intent(in), value :: denominator
    end subroutine adj_vec_pointwisedivide_proc

    subroutine adj_vec_norm_proc(x, norm) bind(c)
      use libadjoint_data_structures
      type(adj_vector), intent(in), value :: x
      adj_scalar_f, intent(out) :: norm
    end subroutine adj_vec_norm_proc

    subroutine adj_vec_dot_product(x, y, val) bind(c)
      use libadjoint_data_structures
      type(adj_vector), intent(in), value :: x, y
      adj_scalar_f, intent(out) :: val
    end subroutine adj_vec_dot_product

    subroutine adj_vec_set_random(x) bind(c)
      use libadjoint_data_structures
      type(adj_vector), intent(inout) :: x
    end subroutine adj_vec_set_random

    subroutine adj_mat_duplicate_proc(matin, matout) bind(c)
      ! Allocate a new matrix, using a given matrix as the model
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_matrix), intent(in), value :: matin
      type(adj_matrix), intent(out) :: matout
    end subroutine adj_mat_duplicate_proc

    subroutine adj_mat_axpy_proc(Y, alpha, X) bind(c)
      ! Computes Y = alpha*X + Y.
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_matrix), intent(inout) :: Y
      adj_scalar_f, intent(in), value :: alpha
      type(adj_matrix), intent(in), value :: X
    end subroutine adj_mat_axpy_proc

    subroutine adj_mat_destroy_proc(mat) bind(c)
      ! Frees space taken by a matrix.
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_matrix), intent(inout) :: mat
    end subroutine adj_mat_destroy_proc

    subroutine adj_nonlinear_colouring_proc(nvar, variables, dependencies, derivative, context, sz, colouring) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      type(adj_variable), intent(in), value :: derivative
      type(c_ptr), intent(in), value :: context
      integer(kind=c_int), intent(in), value :: sz
      integer(kind=c_int), dimension(sz), intent(out) :: colouring
    end subroutine adj_nonlinear_colouring_proc

    subroutine adj_nonlinear_action_proc(nvar, variables, dependencies, input, context, output) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output
    end subroutine adj_nonlinear_action_proc

    subroutine adj_nonlinear_derivative_action_proc(nvar, variables, dependencies, derivative, contraction, hermitian, &
                                                  & input, coefficient, context, output) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      type(adj_variable), intent(in), value :: derivative
      type(adj_vector), intent(in), value :: contraction
      integer(kind=c_int), intent(in), value :: hermitian
      type(adj_vector), intent(in), value :: input
      adj_scalar_f, intent(in), value :: coefficient
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output
    end subroutine adj_nonlinear_derivative_action_proc

    subroutine adj_nonlinear_derivative_assembly_proc(nvar, variables, dependencies, derivative, contraction, hermitian, &
                                                    & context, output) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      integer(kind=c_int), intent(in) :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      type(adj_variable), intent(in) :: derivative
      type(adj_vector), intent(in) :: contraction
      logical(kind=c_bool), intent(in) :: hermitian
      type(c_ptr), intent(in) :: context
      type(adj_matrix), intent(out) :: output
    end subroutine adj_nonlinear_derivative_assembly_proc
   
    subroutine adj_block_action_proc(nvar, variables, dependencies, hermitian, coefficient, input, context, output) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(adj_vector), intent(in), value :: input
      type(c_ptr), intent(in), value :: context
      type(adj_vector), intent(out) :: output
    end subroutine adj_block_action_proc

    subroutine adj_block_assembly_proc(nvar, variables, dependencies, hermitian, coefficient, context, output, rhs) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      integer(kind=c_int), intent(in), value :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      integer(kind=c_int), intent(in), value :: hermitian
      adj_scalar_f, intent(in), value :: coefficient
      type(c_ptr), intent(in), value :: context
      type(adj_matrix), intent(out) :: output
      type(adj_vector), intent(out) :: rhs
    end subroutine adj_block_assembly_proc
    
    subroutine adj_functional_proc(adjointer, timestep, ndepends, dependencies, values, name, output) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_adjointer), intent(in) :: adjointer
      integer(kind=c_int), intent(in), value :: timestep
      integer(kind=c_int), intent(in), value :: ndepends
      type(adj_variable), dimension(ndepends), intent(in) :: dependencies
      type(adj_vector), dimension(ndepends), intent(in) :: values
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
      adj_scalar_f, intent(out) :: output
    end subroutine adj_functional_proc

    subroutine adj_functional_derivative_proc(adjointer, variable, ndepends, dependencies, values, name, output) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_adjointer), intent(in) :: adjointer
      type(adj_variable), intent(in), value :: variable
      integer(kind=c_int), intent(in), value :: ndepends
      type(adj_variable), dimension(ndepends), intent(in) :: dependencies
      type(adj_vector), dimension(ndepends), intent(in) :: values
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
      type(adj_vector), intent(out) :: output
    end subroutine adj_functional_derivative_proc

    subroutine adj_forward_source_proc(adjointer, variable, ndepends, dependencies, values, output, has_output) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_adjointer), intent(in) :: adjointer
      type(adj_variable), intent(in), value :: variable
      integer(kind=c_int), intent(in), value :: ndepends
      type(adj_variable), dimension(ndepends), intent(in) :: dependencies
      type(adj_vector), dimension(ndepends), intent(in) :: values
      type(adj_vector), intent(out) :: output
      integer(kind=c_int), intent(out) :: has_output
    end subroutine adj_forward_source_proc
  end interface

  interface
    function adj_create_variable_c(name, timestep, iteration, auxiliary, var) result(ierr) bind(c, name='adj_create_variable')
      use libadjoint_data_structures
      use iso_c_binding
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
      integer(kind=c_int), intent(in), value :: timestep, iteration, auxiliary
      type(adj_variable), intent(out) :: var
      integer(kind=c_int) :: ierr
    end function adj_create_variable_c

    function adj_variable_get_timestep(var, timestep) result(ierr) bind(c, name='adj_variable_get_timestep')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_variable), intent(in), value :: var
      integer(kind=c_int), intent(out) :: timestep
      integer(kind=c_int) :: ierr
    end function adj_variable_get_timestep

    function adj_variable_get_iteration(var, iteration) result(ierr) bind(c, name='adj_variable_get_iteration')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_variable), intent(in), value :: var
      integer(kind=c_int), intent(out) :: iteration
      integer(kind=c_int) :: ierr
    end function adj_variable_get_iteration

    function adj_variable_set_auxiliary_c(var, auxiliary) result(ierr) bind(c, name='adj_variable_set_auxiliary')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_variable), intent(inout) :: var
      integer(kind=c_int), intent(in), value :: auxiliary
      integer(kind=c_int) :: ierr
    end function adj_variable_set_auxiliary_c

    subroutine adj_chkierr_private_c(ierr, filename, line) bind(c, name='adj_chkierr_private')
      use iso_c_binding
      integer(kind=c_int), intent(in), value :: ierr
      character(kind=c_char), dimension(*), intent(in) :: filename
      integer(kind=c_int), intent(in), value :: line
    end subroutine adj_chkierr_private_c

    function adj_create_nonlinear_block_c(name, ndepends, depends, context, nblock) result(ierr) &
      & bind(c, name='adj_create_nonlinear_block')
      use libadjoint_data_structures
      use iso_c_binding
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
      integer(kind=c_int), intent(in), value :: ndepends
      type(adj_variable), intent(in), dimension(ndepends) :: depends
      type(c_ptr), intent(in), value :: context
      type(adj_nonlinear_block), intent(out) :: nblock
      integer(kind=c_int) :: ierr
    end function adj_create_nonlinear_block_c

    function adj_destroy_nonlinear_block(nblock) result(ierr) bind(c, name='adj_destroy_nonlinear_block')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_nonlinear_block), intent(inout) :: nblock
      integer(kind=c_int) :: ierr
    end function adj_destroy_nonlinear_block

    function adj_nonlinear_block_set_coefficient(nblock, coefficient) result(ierr) &
      & bind(c, name='adj_nonlinear_block_set_coefficient')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_nonlinear_block), intent(inout) :: nblock
      adj_scalar_f, intent(in), value :: coefficient
      integer(kind=c_int) :: ierr
    end function adj_nonlinear_block_set_coefficient

    function adj_create_block_c(name, nblock, context, block) result(ierr) bind(c, name='adj_create_block')
      use libadjoint_data_structures
      use iso_c_binding
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
      ! I want it to be
      ! type(adj_nonlinear_block), intent(in), optional :: nblock
      ! but I can't because optional and bind(c) are incompatible.
      ! So ..
      type(c_ptr), intent(in), value :: nblock
      type(c_ptr), intent(in), value :: context
      type(adj_block), intent(inout) :: block
      integer(kind=c_int) :: ierr
    end function adj_create_block_c

    function adj_destroy_block(block) result(ierr) bind(c, name='adj_destroy_block')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_block), intent(inout) :: block
      integer(kind=c_int) :: ierr
    end function adj_destroy_block

    function adj_block_set_coefficient(block, coefficient) result(ierr) bind(c, name='adj_block_set_coefficient')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_block), intent(inout) :: block
      adj_scalar_f, intent(in), value :: coefficient
      integer(kind=c_int) :: ierr
    end function adj_block_set_coefficient

    function adj_block_set_hermitian_c(block, hermitian) result(ierr) bind(c, name='adj_block_set_hermitian')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_block), intent(inout) :: block
      integer(kind=c_int), intent(in), value :: hermitian
      integer(kind=c_int) :: ierr
    end function adj_block_set_hermitian_c

    function adj_block_set_test_hermitian_c(block, test_hermitian, number_of_tests, tolerance) result(ierr) &
                                        & bind(c, name='adj_block_set_test_hermitian')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_block), intent(inout) :: block
      integer(kind=c_int), intent(in), value :: test_hermitian
      integer(kind=c_int), intent(in), value :: number_of_tests
      adj_scalar_f, intent(in), value :: tolerance
      integer(kind=c_int) :: ierr
    end function adj_block_set_test_hermitian_c

    function adj_nonlinear_block_set_test_hermitian_c(nblock, test_hermitian, number_of_tests, tolerance) result(ierr) &
                                        & bind(c, name='adj_nonlinear_block_set_test_hermitian')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_nonlinear_block), intent(inout) :: nblock
      integer(kind=c_int), intent(in), value :: test_hermitian
      integer(kind=c_int), intent(in), value :: number_of_tests
      adj_scalar_f, intent(in), value :: tolerance
      integer(kind=c_int) :: ierr
    end function adj_nonlinear_block_set_test_hermitian_c
    
    function adj_nonlinear_block_set_test_derivative_c(nblock, test_derivative, number_of_rounds) result(ierr) &
                                        & bind(c, name='adj_nonlinear_block_set_test_derivative')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_nonlinear_block), intent(inout) :: nblock
      integer(kind=c_int), intent(in), value :: test_derivative
      integer(kind=c_int), intent(in), value :: number_of_rounds
      integer(kind=c_int) :: ierr
    end function adj_nonlinear_block_set_test_derivative_c

    function adj_create_equation_c(variable, nblocks, blocks, targets, equation) result(ierr) bind(c, name='adj_create_equation')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_variable), intent(in), value :: variable
      integer(kind=c_int), intent(in), value :: nblocks
      type(adj_block), dimension(nblocks), intent(in) :: blocks
      type(adj_variable), dimension(nblocks), intent(in) :: targets
      type(adj_equation), intent(inout) :: equation
      integer(kind=c_int) :: ierr
    end function adj_create_equation_c

    function adj_equation_set_rhs_dependencies_c(equation, nrhsdeps, rhsdeps, context) result(ierr) &
                                              & bind(c, name='adj_equation_set_rhs_dependencies')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_equation), intent(inout) :: equation
      integer(kind=c_int), intent(in), value :: nrhsdeps
      type(adj_variable), dimension(nrhsdeps), intent(in) :: rhsdeps
      type(c_ptr), intent(in), value :: context
      integer(kind=c_int) :: ierr
    end function adj_equation_set_rhs_dependencies_c

    function adj_destroy_equation(equation) result(ierr) bind(c, name='adj_destroy_equation')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_equation), intent(inout) :: equation
      integer(kind=c_int) :: ierr
    end function adj_destroy_equation

    function adj_create_adjointer(adjointer) result(ierr) bind(c, name='adj_create_adjointer')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int) :: ierr
    end function adj_create_adjointer

    function adj_destroy_adjointer(adjointer) result(ierr) bind(c, name='adj_destroy_adjointer')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int) :: ierr
    end function adj_destroy_adjointer

    function adj_set_option(adjointer, option, choice) result(ierr) bind(c, name='adj_set_option')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: option, choice 
      integer(kind=c_int) :: ierr
    end function adj_set_option

    function adj_equation_count(adjointer, count) result(ierr) bind(c, name='adj_equation_count')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      integer(kind=c_int), intent(inout) :: count
      integer(kind=c_int) :: ierr
    end function adj_equation_count

    function adj_register_equation(adjointer, equation) result(ierr) bind(c, name='adj_register_equation')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      type(adj_equation), intent(in), value :: equation
      integer(kind=c_int) :: ierr
    end function adj_register_equation

    function adj_record_variable(adjointer, variable, storage) result(ierr) bind(c, name='adj_record_variable')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      type(adj_variable), intent(in), value :: variable
      type(adj_storage_data), intent(in), value :: storage
      integer(kind=c_int) :: ierr
    end function adj_record_variable

    function adj_register_operator_callback_c(adjointer, type, name, fnptr) result(ierr) &
           & bind(c, name='adj_register_operator_callback')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: type
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
      type(c_funptr), intent(in), value :: fnptr
      integer(kind=c_int) :: ierr
    end function adj_register_operator_callback_c

    function adj_register_data_callback(adjointer, type, fnptr) result(ierr) bind(c, name='adj_register_data_callback')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: type
      type(c_funptr), intent(in), value :: fnptr
      integer(kind=c_int) :: ierr
    end function adj_register_data_callback
    
    function adj_register_functional_callback_c(adjointer, name, fnptr) &
                                              & result(ierr) bind(c, name='adj_register_functional_callback')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
      type(c_funptr), intent(in), value :: fnptr
      integer(kind=c_int) :: ierr
    end function adj_register_functional_callback_c

    function adj_register_functional_derivative_callback_c(adjointer, name, fnptr) &
                                                      & result(ierr) bind(c, name='adj_register_functional_derivative_callback')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: name
      type(c_funptr), intent(in), value :: fnptr
      integer(kind=c_int) :: ierr
    end function adj_register_functional_derivative_callback_c

    function adj_register_forward_source_callback(adjointer, fnptr) &
                                                      & result(ierr) bind(c, name='adj_register_forward_source_callback')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      type(c_funptr), intent(in), value :: fnptr
      integer(kind=c_int) :: ierr
    end function adj_register_forward_source_callback

    function adj_forget_adjoint_equation(adjointer, equation) result(ierr) bind(c, name='adj_forget_adjoint_equation')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: equation
      integer(kind=c_int) :: ierr
    end function adj_forget_adjoint_equation

    function adj_timestep_count(adjointer, count) result(ierr) bind(c, name='adj_timestep_count')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      integer(kind=c_int), intent(out) :: count
      integer(kind=c_int) :: ierr
    end function adj_timestep_count

    function adj_iteration_count(adjointer, variable, count) result(ierr) bind(c, name='adj_iteration_count')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      type(adj_variable), intent(in), value :: variable
      integer(kind=c_int), intent(out) :: count
      integer(kind=c_int) :: ierr
    end function adj_iteration_count

    function adj_timestep_start_equation(adjointer, timestep, start) result(ierr) bind(c, name='adj_timestep_start_equation')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      integer(kind=c_int), intent(in), value :: timestep
      integer(kind=c_int), intent(out) :: start
      integer(kind=c_int) :: ierr
    end function adj_timestep_start_equation

    function adj_timestep_end_equation(adjointer, timestep, end) result(ierr) bind(c, name='adj_timestep_end_equation')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      integer(kind=c_int), intent(in), value :: timestep
      integer(kind=c_int), intent(out) :: end
      integer(kind=c_int) :: ierr
    end function adj_timestep_end_equation

    function adj_adjointer_check_consistency(adjointer) result(ierr) bind(c, name='adj_adjointer_check_consistency')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      integer(kind=c_int) :: ierr
    end function adj_adjointer_check_consistency

    function adj_timestep_set_times(adjointer, timestep, start, end) result(ierr) bind(c, name='adj_timestep_set_times')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: timestep
      adj_scalar_f, intent(in), value :: start
      adj_scalar_f, intent(in), value :: end
      integer(kind=c_int) :: ierr
    end function adj_timestep_set_times

    function adj_timestep_get_times(adjointer, timestep, start, end) result(ierr) bind(c, name='adj_timestep_get_times')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      integer(kind=c_int), intent(in), value :: timestep
      adj_scalar_f, intent(out) :: start
      adj_scalar_f, intent(out) :: end
      integer(kind=c_int) :: ierr
    end function adj_timestep_get_times

    function adj_timestep_set_functional_dependencies_c(adjointer, timestep, functional, ndepends, dependencies) result(ierr) &
                                                      & bind(c, name='adj_timestep_set_functional_dependencies')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: timestep
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: functional
      integer(kind=c_int), intent(in), value :: ndepends
      type(adj_variable), dimension(*), intent(in) :: dependencies
      integer(kind=c_int) :: ierr
    end function adj_timestep_set_functional_dependencies_c

    function adj_variable_get_ndepending_timesteps_c(adjointer, variable, functional, ntimesteps) result(ierr) &
                                                      & bind(c, name='adj_variable_get_ndepending_timesteps')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      type(adj_variable), intent(in), value :: variable
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: functional
      integer(kind=c_int), intent(out) :: ntimesteps
      integer(kind=c_int) :: ierr
    end function adj_variable_get_ndepending_timesteps_c

    function adj_variable_get_depending_timestep_c(adjointer, variable, functional, i, timesteps) result(ierr) &
                                                      & bind(c, name='adj_variable_get_depending_timestep')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      type(adj_variable), intent(in), value :: variable
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: functional
      integer(kind=c_int), intent(in), value :: i
      integer(kind=c_int), intent(out) :: timesteps
      integer(kind=c_int) :: ierr
    end function adj_variable_get_depending_timestep_c

    function adj_storage_memory_copy(val, mem) result(ierr) bind(c, name='adj_storage_memory_copy')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_vector), intent(in), value :: val
      type(adj_storage_data), intent(inout) :: mem
      integer(kind=c_int) :: ierr
    end function adj_storage_memory_copy
    
    function adj_storage_memory_incref(val, mem) result(ierr) bind(c, name='adj_storage_memory_incref')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_vector), intent(in), value :: val
      type(adj_storage_data), intent(inout) :: mem
      integer(kind=c_int) :: ierr
    end function adj_storage_memory_incref

    function adj_storage_set_compare_c(mem, compare, comparison_tolerance) result(ierr) bind(c, name='adj_storage_set_compare')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_storage_data), intent(inout) :: mem
      integer(kind=c_int), intent(in), value :: compare
      adj_scalar_f, intent(in), value :: comparison_tolerance
      integer(kind=c_int) :: ierr
    end function adj_storage_set_compare_c

    function adj_storage_set_overwrite_c(mem, overwrite) result(ierr) bind(c, name='adj_storage_set_overwrite')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_storage_data), intent(inout) :: mem
      integer(kind=c_int), intent(in), value :: overwrite
      integer(kind=c_int) :: ierr
    end function adj_storage_set_overwrite_c

    function adj_get_adjoint_equation_c(adjointer, equation, functional, lhs, rhs, variable) result(ierr) &
            & bind(c, name='adj_get_adjoint_equation')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: equation
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: functional
      type(adj_matrix), intent(out) :: lhs
      type(adj_vector), intent(out) :: rhs
      type(adj_variable), intent(out) :: variable
      integer(kind=c_int) :: ierr
    end function adj_get_adjoint_equation_c
    
    function adj_get_forward_equation(adjointer, equation, lhs, rhs, variable) result(ierr) &
            & bind(c, name='adj_get_forward_equation')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: equation
      type(adj_matrix), intent(out) :: lhs
      type(adj_vector), intent(out) :: rhs
      type(adj_variable), intent(out) :: variable
      integer(kind=c_int) :: ierr
    end function adj_get_forward_equation
    
    function adj_adjointer_to_html_c(adjointer, filename, type) result(ierr) &
            & bind(c, name='adj_adjointer_to_html')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(in) :: adjointer
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: filename
      integer(kind=c_int), intent(in), value :: type
      integer(kind=c_int) :: ierr
    end function adj_adjointer_to_html_c

    function adj_dict_init(dict) result(ierr) bind(c, name='adj_dict_init')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_dictionary), intent(inout) :: dict
      integer(kind=c_int) :: ierr
    end function adj_dict_init

    function adj_dict_set_c(dict, key, value) result(ierr) bind(c, name='adj_dict_set')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_dictionary), intent(inout) :: dict
      character(kind=c_char), dimension(ADJ_DICT_LEN), intent(in) :: key
      character(kind=c_char), dimension(ADJ_DICT_LEN), intent(in) :: value
      integer(kind=c_int) :: ierr
    end function adj_dict_set_c

    function adj_dict_find_c(dict, key, value) result(ierr) bind(c, name='adj_dict_find')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_dictionary), intent(in) :: dict
      character(kind=c_char), dimension(ADJ_DICT_LEN), intent(in) :: key
      type(c_ptr), intent(out) :: value
      integer(kind=c_int) :: ierr
    end function adj_dict_find_c

    function adj_dict_destroy(dict) result(ierr) bind(c, name='adj_dict_destroy')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_dictionary), intent(inout) :: dict
      integer(kind=c_int) :: ierr
    end function adj_dict_destroy

    subroutine adj_dict_print(dict) bind(c, name='adj_dict_print')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_dictionary), intent(inout) :: dict
    end subroutine adj_dict_print

    function adj_set_petsc_data_callbacks(adjointer) result(ierr) bind(c, name='adj_set_petsc_data_callbacks')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int) :: ierr
    end function adj_set_petsc_data_callbacks

    function adj_evaluate_functional_c(adjointer, timestep, functional, output) result(ierr) &
           & bind(c, name='adj_evaluate_functional')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: timestep
      character(kind=c_char), dimension(ADJ_NAME_LEN), intent(in) :: functional
      adj_scalar_f, intent(out) :: output 
      integer(kind=c_int) :: ierr
    end function adj_evaluate_functional_c

    pure function adj_variable_equal_c(varsa, varsb, sz) result(eq) bind(c, name='adj_variable_equal')
      use libadjoint_data_structures
      use iso_c_binding
      integer(kind=c_int), intent(in), value :: sz
      type(adj_variable), intent(in), dimension(sz) :: varsa, varsb
      integer(kind=c_int) :: eq
    end function adj_variable_equal_c

  end interface

  private :: adj_variable_equal, adj_variable_equal_c
  interface operator (==)
    module procedure adj_variable_equal
  end interface

  contains

  function adj_get_adjoint_equation(adjointer, equation, functional, lhs, rhs, variable) result(ierr)
    type(adj_adjointer), intent(inout) :: adjointer
    integer(kind=c_int), intent(in), value :: equation
    character(len=*), intent(in) :: functional
    type(adj_matrix), intent(out) :: lhs
    type(adj_vector), intent(out) :: rhs
    type(adj_variable), intent(out) :: variable
    integer(kind=c_int) :: ierr
    
    character(kind=c_char), dimension(ADJ_NAME_LEN) :: functional_c
    integer :: j

    do j=1,len_trim(functional)
      functional_c(j) = functional(j:j)
    end do
    do j=len_trim(functional)+1,ADJ_NAME_LEN
      functional_c(j) = c_null_char
    end do
    functional_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_get_adjoint_equation_c(adjointer, equation, functional_c, lhs, rhs, variable)
  end function adj_get_adjoint_equation

  function adj_create_variable(name, timestep, iteration, auxiliary, variable) result(ierr)
    character(len=*), intent(in) :: name
    integer, intent(in) :: timestep, iteration
    logical, intent(in) :: auxiliary
    type(adj_variable), intent(out) :: variable
    integer :: ierr

    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name_c
    integer :: j
    integer :: auxiliary_c

    if (auxiliary) then
      auxiliary_c = ADJ_TRUE
    else
      auxiliary_c = ADJ_FALSE
    end if
    do j=1,len_trim(name)
      name_c(j) = name(j:j)
    end do
    do j=len_trim(name)+1,ADJ_NAME_LEN
      name_c(j) = c_null_char
    end do
    name_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_create_variable_c(name_c, timestep, iteration, auxiliary_c, variable)
  end function adj_create_variable

  function adj_variable_get_name(variable, name) result(ierr)
    type(adj_variable), intent(in) :: variable
    character(len=*), intent(out) :: name
    integer :: ierr
    integer :: j

    name = ' '
    do j=1,min(len(name), ADJ_NAME_LEN)
      if (variable%name(j) /= c_null_char) then
        name(j:j) = variable%name(j)
      else
        name(j:) = ' '
      end if
    end do

    ierr = ADJ_OK
  end function adj_variable_get_name

  function adj_create_nonlinear_block(name, depends, coefficient, context, nblock) result(ierr)
    character(len=*), intent(in) :: name
    type(adj_variable), intent(in), dimension(:) :: depends
    adj_scalar_f, intent(in), optional :: coefficient
    type(c_ptr), intent(in), optional :: context
    type(adj_nonlinear_block), intent(out) :: nblock
    integer :: ierr

    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name_c
    integer :: j
    type(c_ptr) :: context_c

    do j=1,len_trim(name)
      name_c(j) = name(j:j)
    end do
    do j=len_trim(name)+1,ADJ_NAME_LEN
      name_c(j) = c_null_char
    end do
    name_c(ADJ_NAME_LEN) = c_null_char

    if (present(context)) then
      context_c = context
    else
      context_c = c_null_ptr
    endif

    ierr = adj_create_nonlinear_block_c(name_c, size(depends), depends, context_c, nblock)
    if (ierr /= ADJ_OK) return;

    if (present(coefficient)) then
      ierr = adj_nonlinear_block_set_coefficient(nblock, coefficient)
    endif

  end function adj_create_nonlinear_block

  function adj_create_block(name, nblock, context, block) result(ierr)
    use libadjoint_data_structures
    use iso_c_binding
    character(len=*), intent(in) :: name
    type(adj_nonlinear_block), intent(in), optional, target :: nblock
    type(c_ptr), intent(in), optional :: context
    type(adj_block), intent(inout) :: block
    integer :: ierr

    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name_c
    integer :: j
    type(c_ptr) :: nblock_ptr
    type(c_ptr) :: context_c

    do j=1,len_trim(name)
      name_c(j) = name(j:j)
    end do
    do j=len_trim(name)+1,ADJ_NAME_LEN
      name_c(j) = c_null_char
    end do
    name_c(ADJ_NAME_LEN) = c_null_char

    if (present(nblock)) then
      nblock_ptr = c_loc(nblock)
    else
      nblock_ptr = c_null_ptr
    end if

    if (present(context)) then
      context_c = context
    else
      context_c = c_null_ptr
    endif

    ierr = adj_create_block_c(name_c, nblock_ptr, context_c, block)
  end function adj_create_block

  function adj_register_operator_callback(adjointer, type, name, fnptr) result(ierr)
    type(adj_adjointer), intent(inout) :: adjointer
    integer(kind=c_int), intent(in) :: type
    character(len=*), intent(in) :: name
    type(c_funptr), intent(in), value :: fnptr
    integer(kind=c_int) :: ierr

    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name_c
    integer :: j

    do j=1,len_trim(name)
      name_c(j) = name(j:j)
    end do
    do j=len_trim(name)+1,ADJ_NAME_LEN
      name_c(j) = c_null_char
    end do
    name_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_register_operator_callback_c(adjointer, type, name_c, fnptr)
  end function adj_register_operator_callback

  function adj_create_equation(variable, blocks, targets, equation) result(ierr)
    use libadjoint_data_structures
    use iso_c_binding
    type(adj_variable), intent(in), value :: variable
    type(adj_block), dimension(:), intent(in) :: blocks
    type(adj_variable), dimension(:), intent(in) :: targets
    type(adj_equation), intent(inout) :: equation
    integer :: ierr

    if (size(blocks) /= size(targets)) then
      ! Can't set the error message from Fortran, I think?
      ierr = ADJ_ERR_INVALID_INPUTS
      return
    end if

    ierr = adj_create_equation_c(variable, size(blocks), blocks, targets, equation)
  end function adj_create_equation

  function adj_timestep_set_functional_dependencies(adjointer, timestep, functional, dependencies) result(ierr)
    use libadjoint_data_structures
    use iso_c_binding
    type(adj_adjointer), intent(inout) :: adjointer
    integer, intent(in) :: timestep
    character(len=*), intent(in) :: functional
    type(adj_variable), dimension(:), intent(in) :: dependencies
    integer :: ierr

    character(kind=c_char), dimension(ADJ_NAME_LEN) :: functional_c
    integer :: j

    do j=1,len_trim(functional)
      functional_c(j) = functional(j:j)
    end do
    do j=len_trim(functional)+1,ADJ_NAME_LEN
      functional_c(j) = c_null_char
    end do
    functional_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_timestep_set_functional_dependencies_c(adjointer, timestep, functional_c, size(dependencies), dependencies)
  end function adj_timestep_set_functional_dependencies

  function adj_variable_get_ndepending_timesteps(adjointer, variable, functional, ntimesteps) result(ierr)
    type(adj_adjointer), intent(in) :: adjointer
    type(adj_variable), intent(in), value :: variable
    character(len=*), intent(in) :: functional
    integer, intent(out) :: ntimesteps
    integer :: ierr

    character(kind=c_char), dimension(ADJ_NAME_LEN) :: functional_c
    integer :: j

    do j=1,len_trim(functional)
      functional_c(j) = functional(j:j)
    end do
    do j=len_trim(functional)+1,ADJ_NAME_LEN
      functional_c(j) = c_null_char
    end do
    functional_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_variable_get_ndepending_timesteps_c(adjointer, variable, functional_c, ntimesteps)
  end function adj_variable_get_ndepending_timesteps

  function adj_variable_get_depending_timestep(adjointer, variable, functional, i, timestep) result(ierr)
    type(adj_adjointer), intent(in) :: adjointer
    type(adj_variable), intent(in), value :: variable
    character(len=*), intent(in) :: functional
    integer, intent(in), value :: i
    integer, intent(out) :: timestep
    integer :: ierr

    character(kind=c_char), dimension(ADJ_NAME_LEN) :: functional_c
    integer :: j

    do j=1,len_trim(functional)
      functional_c(j) = functional(j:j)
    end do
    do j=len_trim(functional)+1,ADJ_NAME_LEN
      functional_c(j) = c_null_char
    end do
    functional_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_variable_get_depending_timestep_c(adjointer, variable, functional_c, i, timestep)
  end function adj_variable_get_depending_timestep
  
  function adj_register_functional_callback(adjointer, name, fnptr) result(ierr)
    type(adj_adjointer), intent(inout) :: adjointer
    character(len=*), intent(in) :: name
    type(c_funptr), intent(in) :: fnptr
    integer :: ierr

    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name_c
    integer :: j

    do j=1,len_trim(name)
      name_c(j) = name(j:j)
    end do
    do j=len_trim(name)+1,ADJ_NAME_LEN
      name_c(j) = c_null_char
    end do
    name_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_register_functional_callback_c(adjointer, name_c, fnptr)
  end function adj_register_functional_callback

  function adj_register_functional_derivative_callback(adjointer, name, fnptr) result(ierr)
    type(adj_adjointer), intent(inout) :: adjointer
    character(len=*), intent(in) :: name
    type(c_funptr), intent(in) :: fnptr
    integer :: ierr

    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name_c
    integer :: j

    do j=1,len_trim(name)
      name_c(j) = name(j:j)
    end do
    do j=len_trim(name)+1,ADJ_NAME_LEN
      name_c(j) = c_null_char
    end do
    name_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_register_functional_derivative_callback_c(adjointer, name_c, fnptr)
  end function adj_register_functional_derivative_callback

  function adj_adjointer_to_html(adjointer, filename, type) result(ierr) 
    type(adj_adjointer), intent(inout) :: adjointer
    character(len=*), intent(in) :: filename
    integer, intent(in) :: type
    integer :: ierr
    
    character(kind=c_char), dimension(ADJ_NAME_LEN) :: filename_c
    integer :: j

    do j=1,len_trim(filename)
      filename_c(j) = filename(j:j)
    end do
    do j=len_trim(filename)+1,ADJ_NAME_LEN
      filename_c(j) = c_null_char
    end do
    filename_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_adjointer_to_html_c(adjointer, filename_c, type)
  end function adj_adjointer_to_html
  
  function adj_evaluate_functional(adjointer, timestep, functional, output) result(ierr) 
    type(adj_adjointer), intent(inout) :: adjointer
    integer, intent(in) :: timestep
    character(len=*), intent(in) :: functional
    adj_scalar_f, intent(out) :: output
    integer :: ierr
    
    character(kind=c_char), dimension(ADJ_NAME_LEN) :: functional_c
    integer :: j

    do j=1,len_trim(functional)
      functional_c(j) = functional(j:j)
    end do
    do j=len_trim(functional)+1,ADJ_NAME_LEN
      functional_c(j) = c_null_char
    end do
    functional_c(ADJ_NAME_LEN) = c_null_char

    ierr = adj_evaluate_functional_c(adjointer, timestep, functional_c, output)
  end function adj_evaluate_functional
  
  function adj_dict_set(dict, key, value) result(ierr)
    type(adj_dictionary), intent(inout) :: dict
    character(len=*), intent(in) :: key, value
    integer :: ierr

    character(kind=c_char), dimension(ADJ_DICT_LEN) :: key_c
    character(kind=c_char), dimension(ADJ_DICT_LEN) :: value_c

    integer :: j

    key_c = c_null_char
    value_c = c_null_char
    do j=1,len(key)
      key_c(j) = key(j:j)
    end do

    do j=1,len(value)
      value_c(j) = value(j:j)
    end do

    ierr = adj_dict_set_c(dict, key_c, value_c)
  end function adj_dict_set

  function adj_dict_find(dict, key, value) result(ierr)
    type(adj_dictionary), intent(inout) :: dict
    character(len=*), intent(in) :: key
    character(len=*), intent(out) :: value
    integer :: ierr

    character(kind=c_char), dimension(ADJ_DICT_LEN) :: key_c
    character(kind=c_char), dimension(:), pointer :: value_c
    type(c_ptr) :: ptr

    integer :: j

    key_c = c_null_char
    do j=1,len(key)
      key_c(j) = key(j:j)
    end do

    ierr = adj_dict_find_c(dict, key_c, ptr)
    value = " "

    if (ierr == ADJ_OK) then
      call c_f_pointer(ptr, value_c, (/ADJ_DICT_LEN/))
      j = 1
      do while(j <= min(ADJ_DICT_LEN, len(value)) .and. value_c(j) /= c_null_char)
        value(j:j) = value_c(j)
        j = j + 1
      end do
    end if
  end function adj_dict_find

  subroutine adj_chkierr_private(ierr, filename, line)
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: filename
    integer, intent(in) :: line
    character(kind=c_char), dimension(len_trim(filename)+1) :: filename_c
    integer :: j

    do j=1,len_trim(filename)
      filename_c(j) = filename(j:j)
    end do
    filename_c(len_trim(filename)+1) = c_null_char

    call adj_chkierr_private_c(ierr, filename_c, line)
  end subroutine adj_chkierr_private

  subroutine adj_test_assert(bool, testdesc)
    logical, intent(in) :: bool
    character(len=*), intent(in) :: testdesc

    if (.not. bool) then
      print *, "  fail: " // testdesc
    else
      print *, "  pass"
    end if
  end subroutine adj_test_assert

  function adj_variable_set_auxiliary(var, auxiliary) result(ierr)
    type(adj_variable), intent(inout) :: var
    logical, intent(in) :: auxiliary
    integer(kind=c_int) :: auxiliary_c
    integer(kind=c_int) :: ierr

    if (auxiliary) then
      auxiliary_c = ADJ_TRUE
    else
      auxiliary_c = ADJ_FALSE
    end if

    ierr = adj_variable_set_auxiliary_c(var, auxiliary_c)
  end function adj_variable_set_auxiliary

  function adj_equation_set_rhs_dependencies(equation, rhsdeps, context) result(ierr)
    type(adj_equation), intent(inout) :: equation
    type(adj_variable), dimension(:), intent(in), optional :: rhsdeps
    type(c_ptr), intent(in), value, optional :: context
    integer(kind=c_int) :: ierr

    type(c_ptr) :: context_c
    type(adj_variable), dimension(0) :: dummy_rhsdeps
    integer :: nrhsdeps

    if (present(context)) then
      context_c = context
    else
      context_c = c_null_ptr
    end if

    if (present(rhsdeps)) then
      ierr = adj_equation_set_rhs_dependencies_c(equation, size(rhsdeps), rhsdeps, context_c)
    else
      ierr = adj_equation_set_rhs_dependencies_c(equation, 0, dummy_rhsdeps, context_c)
    end if
  end function adj_equation_set_rhs_dependencies

  function adj_storage_set_compare(mem, compare, comparison_tolerance) result(ierr)
    type(adj_storage_data), intent(inout) :: mem
    logical, intent(in) :: compare
    adj_scalar_f, intent(in) :: comparison_tolerance
    integer(kind=c_int) :: ierr

    integer(kind=c_int) :: compare_c

    if (compare) then
      compare_c = ADJ_TRUE
    else
      compare_c = ADJ_FALSE
    end if

    ierr = adj_storage_set_compare_c(mem, compare_c, comparison_tolerance)
  end function adj_storage_set_compare

  function adj_storage_set_overwrite(mem, overwrite) result(ierr)
    type(adj_storage_data), intent(inout) :: mem
    logical, intent(in) :: overwrite
    integer(kind=c_int) :: ierr

    integer(kind=c_int) :: overwrite_c

    if (overwrite) then
      overwrite_c = ADJ_TRUE
    else
      overwrite_c = ADJ_FALSE
    end if

    ierr = adj_storage_set_overwrite_c(mem, overwrite_c)
  end function adj_storage_set_overwrite

  function adj_block_set_test_hermitian(block, test_hermitian, number_of_tests, tolerance) result(ierr)
    use libadjoint_data_structures
    use iso_c_binding
    type(adj_block), intent(inout) :: block
    logical, intent(in) :: test_hermitian
    integer(kind=c_int), intent(in), value :: number_of_tests
    adj_scalar_f, intent(in), value :: tolerance
    integer(kind=c_int) :: ierr

    integer(kind=c_int) :: test_hermitian_c

    if (test_hermitian) then
      test_hermitian_c = ADJ_TRUE
    else
      test_hermitian_c = ADJ_FALSE
    end if

    ierr = adj_block_set_test_hermitian_c(block, test_hermitian_c, number_of_tests, tolerance)
  end function adj_block_set_test_hermitian

  function adj_nonlinear_block_set_test_hermitian(nblock, test_hermitian, number_of_tests, tolerance) result(ierr)
    use libadjoint_data_structures
    use iso_c_binding
    type(adj_nonlinear_block), intent(inout) :: nblock
    logical, intent(in) :: test_hermitian
    integer(kind=c_int), intent(in), value :: number_of_tests
    adj_scalar_f, intent(in), value :: tolerance
    integer(kind=c_int) :: ierr

    integer(kind=c_int) :: test_hermitian_c

    if (test_hermitian) then
      test_hermitian_c = ADJ_TRUE
    else
      test_hermitian_c = ADJ_FALSE
    end if

    ierr = adj_nonlinear_block_set_test_hermitian_c(nblock, test_hermitian_c, number_of_tests, tolerance)
  end function adj_nonlinear_block_set_test_hermitian

  function adj_nonlinear_block_set_test_derivative(nblock, test_derivative, number_of_rounds) result(ierr)
    use libadjoint_data_structures
    use iso_c_binding
    type(adj_nonlinear_block), intent(inout) :: nblock
    logical, intent(in) :: test_derivative
    integer, intent(in) :: number_of_rounds
    integer(kind=c_int) :: ierr

    integer(kind=c_int) :: test_derivative_c

    if (test_derivative) then
      test_derivative_c = ADJ_TRUE
    else
      test_derivative_c = ADJ_FALSE
    end if

    ierr = adj_nonlinear_block_set_test_derivative_c(nblock, test_derivative_c, number_of_rounds)
  end function adj_nonlinear_block_set_test_derivative

  function adj_block_set_hermitian(block, hermitian) result(ierr)
    use libadjoint_data_structures
    use iso_c_binding
    type(adj_block), intent(inout) :: block
    logical, intent(in) :: hermitian
    integer(kind=c_int) :: ierr

    integer(kind=c_int) :: hermitian_c

    if (hermitian) then
      hermitian_c = ADJ_TRUE
    else
      hermitian_c = ADJ_FALSE
    end if
    ierr = adj_block_set_hermitian_c(block, hermitian_c)
  end function adj_block_set_hermitian

  pure function adj_variable_equal(var1, var2)
    type(adj_variable), intent(in) :: var1, var2
    logical :: adj_variable_equal
    type(adj_variable), dimension(1) :: var1tmp, var2tmp

    integer(kind=c_int) :: eq

    var1tmp(1) = var1
    var2tmp(1) = var2
    eq = adj_variable_equal_c(var1tmp, var2tmp, 1)
    adj_variable_equal = (eq == 1)
  end function adj_variable_equal

end module libadjoint

module libadjoint_petsc_data_structures

  use libadjoint
  use iso_c_binding
  implicit none
#ifdef HAVE_PETSC
#include "libadjoint/adj_petsc_f.h"
#endif

  private
#ifdef HAVE_PETSC
  public :: petsc_vec_from_adj_vector
  public :: petsc_vec_to_adj_vector
  public :: petsc_mat_from_adj_matrix
  public :: petsc_mat_to_adj_matrix
  public :: petsc_vec_destroy_proc
  public :: petsc_mat_destroy_proc

  interface petsc_vec_from_adj_vector
    module procedure petsc_vec_from_adj_vector_f
  end interface petsc_vec_from_adj_vector

  interface petsc_vec_to_adj_vector
    module procedure petsc_vec_to_adj_vector_f
  end interface petsc_vec_to_adj_vector

  interface petsc_mat_from_adj_matrix
    module procedure petsc_mat_from_adj_matrix_f
  end interface petsc_mat_from_adj_matrix

  interface petsc_mat_to_adj_matrix
    module procedure petsc_mat_to_adj_matrix_f
  end interface petsc_mat_to_adj_matrix

  interface petsc_vec_destroy_proc
    module procedure petsc_vec_destroy_proc_f
  end interface

  interface petsc_mat_destroy_proc
    module procedure petsc_mat_destroy_proc_f
  end interface
#endif

  contains

  subroutine petsc_vec_destroy_proc_f(x) bind(c)
    ! Destroys a vector.
    use iso_c_binding
    type(adj_vector), intent(inout) :: x
    integer :: ierr
#ifdef HAVE_PETSC
    Vec, pointer :: xvec
    call c_f_pointer(x%ptr, xvec)
    call VecDestroy(xvec, ierr)
    deallocate(xvec)
#endif
  end subroutine petsc_vec_destroy_proc_f

  subroutine petsc_mat_destroy_proc_f(mat) bind(c)
    ! Frees space taken by a matrix.
    use iso_c_binding
    type(adj_matrix), intent(inout) :: mat
    integer :: ierr
#ifdef HAVE_PETSC
    call MatDestroy(petsc_mat_from_adj_matrix(mat), ierr)
#endif
  end subroutine petsc_mat_destroy_proc_f

#ifdef HAVE_PETSC
  function petsc_vec_from_adj_vector_f(input) result(output)
    type(adj_vector), intent(in) :: input
    Vec :: output
    Vec, pointer :: tmp

    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end function petsc_vec_from_adj_vector_f

  function petsc_vec_to_adj_vector_f(input) result(output)
    Vec, intent(in), target :: input
    Vec, pointer :: input_ptr
    type(adj_vector) :: output

    output%klass = 0
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function petsc_vec_to_adj_vector_f

  function petsc_mat_from_adj_matrix_f(input) result(output)
    type(adj_matrix), intent(in) :: input
    Mat :: output
    Mat, pointer :: tmp

    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end function petsc_mat_from_adj_matrix_f

  function petsc_mat_to_adj_matrix_f(input) result(output)
    Mat, intent(in), target :: input
    Mat, pointer :: input_ptr
    type(adj_matrix) :: output

    output%klass = 0
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function petsc_mat_to_adj_matrix_f
#endif
end module libadjoint_petsc_data_structures

