#include "libadjoint/adj_petsc_data_structures.h"

void adj_set_petsc_data_callbacks(adj_adjointer* adjointer)
{
  PetscInitializeNoArguments();
  adj_register_data_callback(adjointer, ADJ_VEC_DUPLICATE_CB, (void (*)(void)) petsc_vec_duplicate_proc);
  adj_register_data_callback(adjointer, ADJ_VEC_AXPY_CB,(void (*)(void)) petsc_vec_axpy_proc);
  adj_register_data_callback(adjointer, ADJ_VEC_DESTROY_CB, (void (*)(void)) petsc_vec_destroy_proc);
  adj_register_data_callback(adjointer, ADJ_MAT_GETVECS_CB, (void (*)(void)) petsc_mat_getvecs_proc);
  adj_register_data_callback(adjointer, ADJ_MAT_AXPY_CB,(void (*)(void)) petsc_mat_axpy_proc);
  adj_register_data_callback(adjointer, ADJ_MAT_DESTROY_CB,(void (*)(void)) petsc_mat_destroy_proc);
  adj_register_data_callback(adjointer, ADJ_MAT_DUPLICATE_CB,(void (*)(void)) petsc_mat_duplicate_proc);
  adj_register_data_callback(adjointer, ADJ_VEC_SETVALUES_CB,(void (*)(void)) petsc_vec_setvalues_proc);
  adj_register_data_callback(adjointer, ADJ_VEC_DIVIDE_CB,(void (*)(void)) petsc_vec_divide_proc);
}

void petsc_vec_duplicate_proc(adj_vector x, adj_vector *newx)
{
#ifdef HAVE_PETSC
  Vec *new_vec=(Vec*) malloc(sizeof(Vec));
  VecDuplicate( *(Vec*) x.ptr, new_vec);
  VecZeroEntries(*new_vec);
  newx->ptr = new_vec;
#endif
}

void petsc_vec_axpy_proc(adj_vector *y, adj_scalar alpha, adj_vector x)
{
#ifdef HAVE_PETSC
    VecAXPY(*(Vec*) y->ptr, alpha, *(Vec*) x.ptr);
#endif
}

void petsc_vec_destroy_proc(adj_vector *x)
{
#ifdef HAVE_PETSC
    VecDestroy(*(Vec*) x->ptr);
#endif
    free(x->ptr);
}

void petsc_mat_duplicate_proc(adj_matrix matin, adj_matrix *matout) 
{
    /* Duplicate a matrix */
#ifdef HAVE_PETSC
    Mat *new_mat=(Mat*) malloc(sizeof(Mat));
    MatDuplicate(*(Mat*) matin.ptr, MAT_DO_NOT_COPY_VALUES, new_mat);
    MatZeroEntries(*new_mat);
    matout->ptr = (adj_matrix *) new_mat;
#endif
}

void petsc_mat_getvecs_proc(adj_matrix mat, adj_vector *left)
{
  /* Get vector(s) compatible with the matrix, i.e. with the same parallel layout */
#ifdef HAVE_PETSC
    MatGetVecs(*(Mat*) mat.ptr, 0, (Vec*) left->ptr);
    VecZeroEntries(*(Vec*) left->ptr);
#endif
}

void petsc_mat_axpy_proc(adj_matrix *Y, adj_scalar alpha, adj_matrix X)
{
    /* Computes Y = alpha*X + Y. */
#ifdef HAVE_PETSC
    MatAXPY(*(Mat*) Y->ptr, alpha, *(Mat*) X.ptr, SAME_NONZERO_PATTERN);
#endif
}

void petsc_mat_destroy_proc(adj_matrix *mat)
{
    /* Frees space taken by a matrix. */
#ifdef HAVE_PETSC
    MatDestroy(*(Mat*) mat->ptr);
#endif
    free(mat->ptr);
}

void petsc_vec_getsize_proc(adj_vector vec, int *sz)
{
#ifdef HAVE_PETSC
  VecGetSize(*(Vec*) vec.ptr, sz);
#endif
}

void petsc_vec_setvalues_proc(adj_vector *vec, adj_scalar scalars[]) 
{
#ifdef HAVE_PETSC
  int i, sz;
  PetscInt* x;
  petsc_vec_getsize_proc(*vec, &sz);
  x = (PetscInt*) malloc(sz * sizeof(PetscInt));
  for (i=0; i<sz;  i++)
    x[i]=i;
  VecSetValues(*(Vec*) vec->ptr, sz, x, scalars, INSERT_VALUES);
  VecAssemblyBegin(*(Vec*) vec->ptr);
  VecAssemblyEnd(*(Vec*) vec->ptr);
  free(x);
#endif
}

void petsc_vec_divide_proc(adj_vector numerator, adj_vector denominator, adj_vector *output)
{
#ifdef HAVE_PETSC
  VecPointwiseDivide(*(Vec*) output->ptr, *(Vec*) numerator.ptr, *(Vec*) denominator.ptr);
#endif
}

#ifdef HAVE_PETSC
adj_vector petsc_vec_to_adj_vector(Vec* v)
{
  adj_vector vv;
  vv.ptr = (void*)v;
  vv.klass = 0;
  return vv;
}

Vec petsc_vec_from_adj_vector(adj_vector vv)
{
  Vec v;
  v = *(Vec*) vv.ptr;
  return v;
}
#endif
