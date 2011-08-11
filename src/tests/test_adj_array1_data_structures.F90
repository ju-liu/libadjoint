#include "libadjoint/adj_fortran.h"

subroutine test_adj_array1_data_structures

  use libadjoint
  use libadjoint_array1_data_structures
  implicit none

  adj_scalar_f, dimension(2) :: rank1_array 
  adj_scalar_f, dimension(2,2) :: rank2_array
  adj_scalar_f, dimension(:, :), pointer :: rank2_array_2
  type(adj_vector) :: vec
  type(adj_matrix) :: mat1, mat2
  integer :: sz
  adj_scalar_f :: val


  ! Get size
  vec = array1_vec_to_adj_vector(rank1_array)
  call array1_vec_get_size(vec, sz)
  call adj_test_assert(sz == 2, "Vector size should be 15")
  
  ! Norm
  rank1_array(:) = 2
  call array1_vec_get_norm(vec, val)
  call adj_test_assert(abs(val-2.8284271247461903)<1e-7, "Vector norm incorrect")
 
  ! Dot product
  call array1_vec_dot_product(vec, vec, val)
  call adj_test_assert(val==8, "Vector dot product incorrect")

  ! Get random
  call array1_vec_set_random(vec) 

  ! Duplicate mat
  mat1 = array1_mat_to_adj_matrix(rank2_array)
  call array1_mat_duplicate(mat1, mat2)
  rank2_array_2 => array1_mat_from_adj_matrix(mat2)

  call adj_test_assert(maxval(abs(rank2_array_2(:, :))) == 0.0, 'New matrix should be zero')

  ! Mat axpy
  rank2_array_2 = 1.0
  val = 2.0
  call array1_mat_axpy(mat2, val, mat2)
  call adj_test_assert(maxval(abs(rank2_array_2(:, :))) == 3.0, 'Wrong error in mat_axpy')

  ! Mat destroy
  call array1_mat_destroy(mat2)

end subroutine test_adj_array1_data_structures
