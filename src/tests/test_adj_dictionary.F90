subroutine test_adj_dictionary
#include "adj_fortran.h"
  use libadjoint
  implicit none

  type(adj_dictionary) :: dict
  character(len=ADJ_DICT_LEN) :: value
  integer :: ierr

  ierr = adj_dict_init(dict)
  call adj_test_assert(ierr == ADJ_OK, "Should always work")

  ierr = adj_dict_set(dict, "Hello", "World")
  call adj_test_assert(ierr == ADJ_OK, "Should always work")
  ierr = adj_dict_find(dict, "Hello", value)
  call adj_test_assert(ierr == ADJ_OK, "Should have found it")
  call adj_test_assert(trim(value) == "World", "Should be World")

  ierr = adj_dict_set(dict, "Goodbye", "Martha")
  ierr = adj_dict_find(dict, "Hello", value)
  call adj_test_assert(ierr == ADJ_OK, "Should have found it")
  call adj_test_assert(trim(value) == "World", "Should be World")
  ierr = adj_dict_find(dict, "Goodbye", value)
  call adj_test_assert(ierr == ADJ_OK, "Should have found it")
  call adj_test_assert(trim(value) == "Martha", "Should be Martha")

  ierr = adj_dict_find(dict, "Adios", value)
  call adj_test_assert(ierr == ADJ_ERR_DICT_FAILED, "No such key")

  ierr = adj_dict_destroy(dict)
  call adj_test_assert(ierr == ADJ_OK, "Should always work")
end subroutine test_adj_dictionary
