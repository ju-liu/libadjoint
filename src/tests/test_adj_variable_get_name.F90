#include "libadjoint/adj_fortran.h"

subroutine test_adj_variable_get_name

  use libadjoint
  implicit none

  integer :: ierr
  character(len=ADJ_NAME_LEN) :: name
  type(adj_variable) :: var

  ierr = adj_create_variable("Velocity", 17, 0, .false., var)
  call adj_test_assert(ierr == 0, "Should have worked")

  ierr = adj_variable_get_name(var, name)
  call adj_test_assert(trim(name) == "Velocity", "Should be Velocity")
end subroutine test_adj_variable_get_name
