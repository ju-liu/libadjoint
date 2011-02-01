#include "libadjoint/adj_fortran.h"

subroutine test_adj_variable_get_timestep

  use libadjoint
  implicit none

  integer :: ierr
  type(adj_variable) :: var
  integer :: timestep

  ierr = adj_create_variable("Velocity", 17, 0, 0, var)
  call adj_test_assert(ierr == 0, "Should have worked")

  ierr = adj_variable_get_timestep(var, timestep)
  call adj_test_assert(timestep == 17, "Should be 17")
end subroutine test_adj_variable_get_timestep
