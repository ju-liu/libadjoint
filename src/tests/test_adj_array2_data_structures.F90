#include "libadjoint/adj_fortran.h"

subroutine test_adj_array2_data_structures

  use libadjoint
  use libadjoint_array2_data_structures
  implicit none

  adj_scalar_f, dimension(4) :: rank1_array
  adj_scalar_f, dimension(2, 2) :: rank2_array, rank2_array_3 
  adj_scalar_f, dimension(:,:), pointer :: rank2_array_2 
  adj_scalar_f, dimension(2,2,2,2) :: rank4_array
  adj_scalar_f, dimension(:, :,:,:), pointer :: rank4_array_2
  type(adj_vector) :: vec, vec2
  type(adj_matrix) :: mat1, mat2
  integer :: sz
  adj_scalar_f :: val


  ! Get size
  rank2_array = 1.0
  vec = array2_vec_to_adj_vector(rank2_array)
  call array2_vec_get_size(vec, sz)
  call adj_test_assert(sz == 4, "Vector size should be 4")
  
  ! Vec duplicate
  call array2_vec_duplicate(vec, vec2)
  rank2_array_2 => array2_vec_from_adj_vector(vec2)
  call adj_test_assert(maxval(abs(rank2_array_2(:,:))) == 0.0, 'New vector should be zero') 
  
  ! Vec axpy
  rank2_array_2 = 5.0
  val = 0.5
  call array2_vec_axpy(vec2, val, vec)
  call adj_test_assert(maxval(abs(rank2_array_2(:,:)-5.5)) == 0.0, 'wrong value after vec axpy')

  ! Vec set values 
  rank1_array = 2.0
  call array2_vec_set_values(vec, rank1_array)

  ! Norm
  call array2_vec_get_norm(vec, val)
  call adj_test_assert(abs(val-4.0)<1e-7, "Vector norm incorrect")
 
  ! Dot product
  call array2_vec_dot_product(vec, vec, val)
  call adj_test_assert(val==16, "Vector dot product incorrect")

  ! Get random
  call array2_vec_set_random(vec) 

  ! Vec destroy
  call array2_vec_destroy(vec2)

  ! Duplicate mat
  mat1 = array2_mat_to_adj_matrix(rank4_array)
  call array2_mat_duplicate(mat1, mat2)
  rank4_array_2 => array2_mat_from_adj_matrix(mat2)

  call adj_test_assert(maxval(abs(rank4_array_2(:,:,:,:))) == 0.0, 'New matrix should be zero')

  ! Mat axpy
  rank4_array_2 = 1.0
  val = 2.0
  call array2_mat_axpy(mat2, val, mat2)
  call adj_test_assert(maxval(abs(rank4_array_2(:,:,:,:))) == 3.0, 'Wrong value after mat_axpy')

  ! Mat destroy
  call array2_mat_destroy(mat2)

end subroutine test_adj_array2_data_structures
