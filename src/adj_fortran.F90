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
    integer(kind=c_int) :: functional
  end type adj_variable

  type, bind(c) :: adj_nonlinear_block
    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name
    adj_scalar_f :: coefficient
    integer(kind=c_int) :: ndepends
    type(c_ptr) :: depends
    type(c_ptr) :: context
  end type adj_nonlinear_block

  type, bind(c) :: adj_block
    character(kind=c_char), dimension(ADJ_NAME_LEN) :: name
    integer(kind=c_int) :: has_nonlinear_block
    type(adj_nonlinear_block) :: nonlinear_block
    type(c_ptr) :: context
    integer(kind=c_int) :: hermitian
    adj_scalar_f :: coefficient
  end type adj_block

  type, bind(c) :: adj_equation
    type(adj_variable) :: variable
    integer(kind=c_int) :: nblocks
    type(c_ptr) :: blocks
    type(c_ptr) :: targets
    integer(kind=c_int) :: nrhsdeps
    type(c_ptr) :: rhsdeps
  end type adj_equation

  type, bind(c) :: adj_data_callbacks
    type(c_funptr) :: vec_duplicate
    type(c_funptr) :: vec_axpy
    type(c_funptr) :: vec_destroy
    type(c_funptr) :: vec_setvalues
    type(c_funptr) :: vec_getsize
    type(c_funptr) :: vec_divide

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
    type(adj_vector) :: value
    type(c_ptr) :: filename
  end type adj_storage_data
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

    subroutine adj_vec_setvalues_proc(vec, scalars) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      type(adj_vector), intent(inout) :: vec
      adj_scalar_f, dimension(*), intent(in) :: scalars
    end subroutine adj_vec_setvalues_proc

    subroutine adj_vec_pointwisedivide_proc(numerator, denominator) bind(c)
      use libadjoint_data_structures
      type(adj_vector), intent(inout) :: numerator
      type(adj_vector), intent(in), value :: denominator
    end subroutine adj_vec_pointwisedivide_proc

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

    subroutine adj_nonlinear_colouring_proc(nvar, variables, dependencies, derivative, sz, colouring) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      integer(kind=c_int), intent(in) :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      type(adj_variable), intent(in) :: derivative
      integer(kind=c_int), intent(in) :: sz
      integer(kind=c_int), dimension(sz), intent(out) :: colouring
    end subroutine adj_nonlinear_colouring_proc

    subroutine adj_nonlinear_action_proc(nvar, variables, dependencies, input, context, output) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      integer(kind=c_int), intent(in) :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      type(adj_vector), intent(in) :: input
      type(c_ptr), intent(in) :: context
      type(adj_vector), intent(out) :: output
    end subroutine adj_nonlinear_action_proc

    subroutine adj_nonlinear_derivative_action_proc(nvar, variables, dependencies, derivative, contraction, hermitian, &
                                                  & input, context, output) bind(c)
      use iso_c_binding
      use libadjoint_data_structures
      integer(kind=c_int), intent(in) :: nvar
      type(adj_variable), dimension(nvar), intent(in) :: variables
      type(adj_vector), dimension(nvar), intent(in) :: dependencies
      type(adj_variable), intent(in) :: derivative
      type(adj_vector), intent(in) :: contraction
      logical(kind=c_bool), intent(in) :: hermitian
      type(adj_vector), intent(in) :: input
      type(c_ptr), intent(in) :: context
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

    function adj_set_rhs_dependencies(equation, nrhsdeps, rhsdeps) result(ierr) bind(c, name='adj_set_rhs_dependencies')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_equation), intent(inout) :: equation
      integer(kind=c_int), intent(in), value :: nrhsdeps
      type(adj_variable), dimension(nrhsdeps), intent(in) :: rhsdeps
      integer(kind=c_int) :: ierr
    end function adj_set_rhs_dependencies

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
      type(adj_adjointer), intent(inout) :: adjointer
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
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(out) :: count
      integer(kind=c_int) :: ierr
    end function adj_timestep_count

    function adj_timestep_start(adjointer, timestep, start) result(ierr) bind(c, name='adj_timestep_start')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: timestep
      integer(kind=c_int), intent(out) :: start
      integer(kind=c_int) :: ierr
    end function adj_timestep_start

    function adj_timestep_end(adjointer, timestep, end) result(ierr) bind(c, name='adj_timestep_end')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: timestep
      integer(kind=c_int), intent(out) :: end
      integer(kind=c_int) :: ierr
    end function adj_timestep_end

    function adj_storage_memory(val) result(mem) bind(c, name='adj_storage_memory')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_vector), intent(in), value :: val
      type(adj_storage_data) :: mem
    end function adj_storage_memory

    function adj_get_adjoint_equation(adjointer, equation, functional, lhs, rhs, variable) result(ierr) &
            & bind(c, name='adj_get_adjoint_equation')
      use libadjoint_data_structures
      use iso_c_binding
      type(adj_adjointer), intent(inout) :: adjointer
      integer(kind=c_int), intent(in), value :: equation
      integer(kind=c_int), intent(in), value :: functional
      type(adj_matrix), intent(out) :: lhs
      type(adj_vector), intent(out) :: rhs
      type(adj_variable), intent(out) :: variable
      integer(kind=c_int) :: ierr
    end function adj_get_adjoint_equation
  end interface

  contains

  function adj_create_variable(name, timestep, iteration, auxiliary, variable) result(ierr)
    character(len=*), intent(in) :: name
    integer, intent(in) :: timestep, iteration
    integer, intent(in) :: auxiliary
    type(adj_variable), intent(out) :: variable
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

    ierr = adj_create_variable_c(name_c, timestep, iteration, auxiliary, variable)
  end function adj_create_variable

  function adj_variable_get_name(variable, name) result(ierr)
    type(adj_variable), intent(in) :: variable
    character(len=ADJ_NAME_LEN), intent(out) :: name
    integer :: ierr
    integer :: j

    do j=1,ADJ_NAME_LEN
      if (variable%name(j) /= c_null_char) then
        name(j:j) = variable%name(j)
      else
        name(j:) = ' '
        exit
      end if
    end do

    ierr = ADJ_ERR_OK
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
    if (ierr /= ADJ_ERR_OK) return;

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

end module libadjoint

module libadjoint_petsc_data_structures

  use libadjoint
  use iso_c_binding
  implicit none
#ifdef HAVE_PETSC
#include "libadjoint/adj_petsc_f.h"
#endif

  private
  public :: adj_set_petsc_data_callbacks
#ifdef HAVE_PETSC
  public :: petsc_vec_from_adj_vector
  public :: petsc_vec_to_adj_vector
  public :: petsc_mat_from_adj_matrix
  public :: petsc_mat_to_adj_matrix
  public :: petsc_vec_destroy_proc
  public :: petsc_mat_destroy_proc
#endif

  contains

  function adj_set_petsc_data_callbacks(adjointer) result(ierr)
    type(adj_adjointer), intent(inout) :: adjointer
    integer(kind=c_int) :: ierr

#ifdef HAVE_PETSC
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
#endif

    ierr = adj_register_data_callback(adjointer, ADJ_VEC_DUPLICATE_CB, c_funloc(petsc_vec_duplicate_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_AXPY_CB, c_funloc(petsc_vec_axpy_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_DESTROY_CB, c_funloc(petsc_vec_destroy_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_SETVALUES_CB, c_funloc(petsc_vec_setvalues_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_VEC_DIVIDE_CB, c_funloc(petsc_vec_divide_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_MAT_AXPY_CB, c_funloc(petsc_mat_axpy_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_MAT_DESTROY_CB, c_funloc(petsc_mat_destroy_proc))
    call adj_chkierr(ierr)
    ierr = adj_register_data_callback(adjointer, ADJ_MAT_DUPLICATE_CB, c_funloc(petsc_mat_duplicate_proc))
    call adj_chkierr(ierr)
  end function adj_set_petsc_data_callbacks

  subroutine petsc_vec_duplicate_proc(x, newx) bind(c)
    type(adj_vector), intent(in), value :: x
    type(adj_vector), intent(out) :: newx
    integer :: ierr
#ifdef HAVE_PETSC
    Vec :: x_petsc
    Vec :: newx_petsc
    x_petsc = petsc_vec_from_adj_vector(x)
    call VecDuplicate(x_petsc, newx_petsc, ierr)
    call VecZeroEntries(newx_petsc, ierr)
    newx = petsc_vec_to_adj_vector(newx_petsc)
#endif
  end subroutine petsc_vec_duplicate_proc

  subroutine petsc_vec_axpy_proc(y, alpha, x) bind(c)
    type(adj_vector), intent(inout) :: y
    adj_scalar_f, intent(in), value :: alpha
    type(adj_vector), intent(in), value :: x
    integer :: ierr
#ifdef HAVE_PETSC
    Vec :: y_petsc, x_petsc
    y_petsc = petsc_vec_from_adj_vector(y)
    x_petsc = petsc_vec_from_adj_vector(x)
    call VecAXPY(y_petsc, alpha, x_petsc, ierr)
#endif
  end subroutine petsc_vec_axpy_proc

  subroutine petsc_vec_destroy_proc(x) bind(c)
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
  end subroutine petsc_vec_destroy_proc

  subroutine petsc_vec_setvalues_proc(vec, scalars) bind(c)
    type(adj_vector), intent(inout) :: vec
    adj_scalar_f, dimension(*), intent(in) :: scalars
    integer :: ierr
    integer :: i
    integer :: sz
#ifdef HAVE_PETSC
    call VecGetSize(petsc_vec_from_adj_vector(vec), sz, ierr)
    call VecSetValues(petsc_vec_from_adj_vector(vec), sz, (/ (i, i=0,sz-1) /), scalars, INSERT_VALUES, ierr)
    call VecAssemblyBegin(petsc_vec_from_adj_vector(vec), ierr)
    call VecAssemblyEnd(petsc_vec_from_adj_vector(vec), ierr)
#endif
  end subroutine petsc_vec_setvalues_proc

  subroutine petsc_vec_divide_proc(numerator, denominator) bind(c)
    ! Computes numerator = numerator/denominator pointwise
    type(adj_vector), intent(inout) :: numerator 
    type(adj_vector), intent(in), value :: denominator
    integer :: ierr
#ifdef HAVE_PETSC

    call VecPointwiseDivide(petsc_vec_from_adj_vector(numerator), petsc_vec_from_adj_vector(numerator), &
                          & petsc_vec_from_adj_vector(denominator), ierr)
#endif
  end subroutine petsc_vec_divide_proc

  subroutine petsc_mat_duplicate_proc(matin, matout) bind(c)
    ! Duplicate a matrix
    use iso_c_binding
    type(adj_matrix), intent(in), value :: matin
    type(adj_matrix), intent(out) :: matout
    integer :: ierr
#ifdef HAVE_PETSC
    Mat :: matout_petsc
    call MatDuplicate(petsc_mat_from_adj_matrix(matin), MAT_DO_NOT_COPY_VALUES, matout_petsc, ierr)
    call MatZeroEntries(matout_petsc, ierr)
    matout = petsc_mat_to_adj_matrix(matout_petsc)
#endif
  end subroutine petsc_mat_duplicate_proc

  subroutine petsc_mat_axpy_proc(Y, alpha, X) bind(c)
    ! Computes Y = alpha*X + Y.
    use iso_c_binding
    type(adj_matrix), intent(inout) :: Y
    adj_scalar_f, intent(in), value :: alpha
    type(adj_matrix), intent(in), value :: X
    integer :: ierr
#ifdef HAVE_PETSC
    call MatAXPY(petsc_mat_from_adj_matrix(Y), alpha, petsc_mat_from_adj_matrix(X), SAME_NONZERO_PATTERN, ierr)
#endif
  end subroutine petsc_mat_axpy_proc

  subroutine petsc_mat_destroy_proc(mat) bind(c)
    ! Frees space taken by a matrix.
    use iso_c_binding
    type(adj_matrix), intent(inout) :: mat
    integer :: ierr
#ifdef HAVE_PETSC
    call MatDestroy(petsc_mat_from_adj_matrix(mat), ierr)
#endif
  end subroutine petsc_mat_destroy_proc

#ifdef HAVE_PETSC
  function petsc_vec_from_adj_vector(input) result(output)
    type(adj_vector), intent(in) :: input
    Vec :: output
    Vec, pointer :: tmp

    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end function petsc_vec_from_adj_vector

  function petsc_vec_to_adj_vector(input) result(output)
    Vec, intent(in), target :: input
    Vec, pointer :: input_ptr
    type(adj_vector) :: output

    output%klass = 0
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function petsc_vec_to_adj_vector

  function petsc_mat_from_adj_matrix(input) result(output)
    type(adj_matrix), intent(in) :: input
    Mat :: output
    Mat, pointer :: tmp

    call c_f_pointer(input%ptr, tmp)
    output = tmp
  end function petsc_mat_from_adj_matrix

  function petsc_mat_to_adj_matrix(input) result(output)
    Mat, intent(in), target :: input
    Mat, pointer :: input_ptr
    type(adj_matrix) :: output

    output%klass = 0
    allocate(input_ptr)
    input_ptr = input
    output%ptr = c_loc(input_ptr)
  end function petsc_mat_to_adj_matrix
#endif
end module libadjoint_petsc_data_structures

