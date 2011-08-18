#include "libadjoint/adj_fortran.h"

subroutine test_adj_datatype_size

  use libadjoint
  use libadjoint_test_tools
  implicit none

  type(adj_adjointer) :: adjointer
  integer :: ierr, size_f
  integer(kind=c_int) :: size_c

  ierr = adj_create_adjointer(adjointer)
  size_f = adj_sizeof_adjointer()
  size_c = sizeof(adjointer)
  if (size_c /= size_f) then
    print *, "Sizeof adjointer in C: ", size_c
    print *, "Sizeof adjointer in Fortran: ", size_f
  end if
  call adj_test_assert(size_f == size_c, "Sizeof adjointer differ between C and Fortran. Check your Fortran interface.")

end subroutine test_adj_datatype_size
