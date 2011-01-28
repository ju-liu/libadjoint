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
end module libadjoint_data_structures

module libadjoint

  use libadjoint_data_structures
  use iso_c_binding
  implicit none


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
      character(kind=c_char), dimension(*), intent(in) :: filename, line
    end subroutine adj_chkierr_private_c
  end interface

  contains

  function adj_create_variable(name, timestep, iteration, auxiliary, var) result(ierr)
    character(len=*), intent(in) :: name
    integer, intent(in) :: timestep, iteration
    integer, intent(in) :: auxiliary
    type(adj_variable), intent(out) :: var
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

    ierr = adj_create_variable_c(name_c, timestep, iteration, auxiliary, var)
  end function adj_create_variable

  function adj_variable_get_name(var, name) result(ierr)
    type(adj_variable), intent(in) :: var
    character(len=ADJ_NAME_LEN), intent(out) :: name
    integer :: ierr
    integer :: j

    do j=1,ADJ_NAME_LEN
      if (var%name(j) /= c_null_char) then
        name(j:j) = var%name(j)
      else
        name(j:) = ' '
        exit
      end if
    end do

    ierr = ADJ_ERR_OK
  end function adj_variable_get_name

  subroutine adj_chkierr_private(ierr, filename, line)
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: filename, line
    character(kind=c_char), dimension(len_trim(filename)+1) :: filename_c
    character(kind=c_char), dimension(len_trim(line)+1) :: line_c
    integer :: j

    do j=1,len_trim(filename)
      filename_c(j) = filename(j:j)
    end do
    filename_c(len_trim(filename)+1) = c_null_char

    do j=1,len_trim(line)
      line_c(j) = line(j:j)
    end do
    line_c(len_trim(line)+1) = c_null_char

    call adj_chkierr_private_c(ierr, filename_c, line_c)
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
