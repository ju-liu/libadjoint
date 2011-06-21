#include "libadjoint/adj_fortran.h"
subroutine test_adj_timestep_count

  use libadjoint
  implicit none

  type(adj_equation) :: equation
  type(adj_adjointer) :: adjointer
  type(adj_variable) :: u
  type(adj_block) :: I
  integer :: ierr
  integer :: count, start, end

  ierr = adj_create_adjointer(adjointer)
  ierr = adj_timestep_count(adjointer, count)
  call adj_test_assert(count == 0, "Should be zero")

  ierr = adj_create_variable("Velocity", timestep=23, iteration=0, auxiliary=.false., variable=u)
  call adj_chkierr(ierr)
  ierr = adj_create_block("Identity", block=I)
  call adj_chkierr(ierr)
  ierr = adj_create_equation(variable=u, blocks=(/I/), targets=(/u/), equation=equation)
  call adj_chkierr(ierr)
  ierr = adj_register_equation(adjointer, equation)
  call adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "Shouldn't work because the first equation must have timestep=0.")
  ierr = adj_destroy_equation(equation)
  call adj_chkierr(ierr)

  ierr = adj_create_variable("Velocity", timestep=0, iteration=0, auxiliary=.false., variable=u)
  call adj_chkierr(ierr)
  ierr = adj_create_equation(variable=u, blocks=(/I/), targets=(/u/), equation=equation)
  call adj_chkierr(ierr)
  ierr = adj_register_equation(adjointer, equation)
  call adj_chkierr(ierr)
  ierr = adj_destroy_equation(equation)
  call adj_chkierr(ierr)
  ierr = adj_create_variable("Velocity", timestep=0, iteration=1, auxiliary=.false., variable=u)
  call adj_chkierr(ierr)
  ierr = adj_create_equation(variable=u, blocks=(/I/), targets=(/u/), equation=equation)
  call adj_chkierr(ierr)
  ierr = adj_register_equation(adjointer, equation)
  call adj_chkierr(ierr)
  ierr = adj_destroy_equation(equation)
  call adj_chkierr(ierr)

  ierr = adj_create_variable("Velocity", timestep=2, iteration=0, auxiliary=.false., variable=u)
  call adj_chkierr(ierr)
  ierr = adj_create_equation(variable=u, blocks=(/I/), targets=(/u/), equation=equation)
  call adj_chkierr(ierr)
  ierr = adj_register_equation(adjointer, equation)
  call adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "Shouldn't work because this equation must have timestep=0 or timestep=1.")
  ierr = adj_destroy_equation(equation)
  call adj_chkierr(ierr)

  ierr = adj_create_variable("Velocity", timestep=1, iteration=0, auxiliary=.false., variable=u)
  call adj_chkierr(ierr)
  ierr = adj_create_equation(variable=u, blocks=(/I/), targets=(/u/), equation=equation)
  call adj_chkierr(ierr)
  ierr = adj_register_equation(adjointer, equation)
  call adj_chkierr(ierr)
  ierr = adj_destroy_equation(equation)
  call adj_chkierr(ierr)

  ierr = adj_timestep_count(adjointer, count)
  call adj_test_assert(count == 2, "Should be 2")
  ierr = adj_timestep_start_equation(adjointer, timestep=2, start=start)
  call adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "No such timestep 2")
  ierr = adj_timestep_end_equation(adjointer, timestep=2, end=end)
  call adj_test_assert(ierr == ADJ_ERR_INVALID_INPUTS, "No such timestep 2")

  ierr = adj_timestep_start_equation(adjointer, timestep=0, start=start)
  call adj_test_assert(start == 0, "Start should be zero")
  ierr = adj_timestep_start_equation(adjointer, timestep=1, start=start)
  call adj_test_assert(start == 2, "Start should be two")

  ierr = adj_timestep_end_equation(adjointer, timestep=0, end=end)
  call adj_test_assert(end == 1, "End should be one")
  ierr = adj_timestep_end_equation(adjointer, timestep=1, end=end)
  call adj_test_assert(end == 2, "End should be two")
end subroutine test_adj_timestep_count
