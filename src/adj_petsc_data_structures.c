#include "libadjoint/adj_petsc_data_structures.h"

int adj_set_petsc_data_callbacks(adj_adjointer* adjointer)
{
  int ierr;
#ifdef HAVE_PETSC
  PetscInitializeNoArguments();
  ierr = adj_register_data_callback(adjointer, ADJ_VEC_DUPLICATE_CB, (void (*)(void)) petsc_vec_duplicate_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_VEC_AXPY_CB,(void (*)(void)) petsc_vec_axpy_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_VEC_DESTROY_CB, (void (*)(void)) petsc_vec_destroy_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_VEC_SET_VALUES_CB,(void (*)(void)) petsc_vec_setvalues_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_VEC_DIVIDE_CB,(void (*)(void)) petsc_vec_divide_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_VEC_GET_NORM_CB,(void (*)(void)) petsc_vec_getnorm_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_VEC_SET_RANDOM_CB,(void (*)(void)) petsc_vec_set_random_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_VEC_DOT_PRODUCT_CB,(void (*)(void)) petsc_vec_dot_product_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_MAT_AXPY_CB,(void (*)(void)) petsc_mat_axpy_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_MAT_DESTROY_CB,(void (*)(void)) petsc_mat_destroy_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_MAT_DUPLICATE_CB,(void (*)(void)) petsc_mat_duplicate_proc);
  adj_chkierr(ierr);
  ierr = adj_register_data_callback(adjointer, ADJ_SOLVE_CB,(void (*)(void)) petsc_solve_proc);
  adj_chkierr(ierr);
#else
  ierr = ADJ_ERR_INVALID_INPUTS;
  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Sorry, libadjoint was compiled without PETSc support.");
#endif

  return ierr;
}

void petsc_vec_duplicate_proc(adj_vector x, adj_vector *newx)
{
#ifdef HAVE_PETSC
  Vec *new_vec=(Vec*) malloc(sizeof(Vec));
  VecDuplicate( *(Vec*) x.ptr, new_vec);
  VecZeroEntries(*new_vec);
  newx->ptr = new_vec;
#else
  (void) x;
  (void) newx;
#endif
}

void petsc_vec_axpy_proc(adj_vector *y, adj_scalar alpha, adj_vector x)
{
#ifdef HAVE_PETSC
    VecAXPY(*(Vec*) y->ptr, alpha, *(Vec*) x.ptr);
#else
    (void) y;
    (void) alpha;
    (void) x;
#endif
}

void petsc_vec_destroy_proc(adj_vector *x)
{
#ifdef HAVE_PETSC
    VecDestroy(*(Vec*) x->ptr);
    free(x->ptr);
#else
    (void) x;
#endif
}

void petsc_vec_getnorm_proc(adj_vector x, adj_scalar* norm)
{
#ifdef HAVE_PETSC
    PetscScalar petsc_norm;
    VecNorm(*(Vec*) x.ptr, NORM_2, &petsc_norm);
    *norm = (adj_scalar) petsc_norm;
#else
    (void) x;
    (void) norm;
#endif
}

void petsc_vec_dot_product_proc(adj_vector x, adj_vector y, adj_scalar* val)
{
#ifdef HAVE_PETSC
  PetscScalar petsc_val;
  VecDot(petsc_vec_from_adj_vector(x), petsc_vec_from_adj_vector(y), &petsc_val);
  *val = (adj_scalar) petsc_val;
#else
    (void) x;
    (void) y;
    (void) val;
#endif
}

void petsc_vec_set_random_proc(adj_vector* x)
{
#ifdef HAVE_PETSC
    struct timeval tval;

    PetscRandom rctx;
    PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
    PetscRandomSetType(rctx, PETSCRAND48);

    gettimeofday(&tval, NULL);
    /* XOR the microseconds since the last whole second since the epoch with the PID of the
       current process. That should make it random enough */
    PetscRandomSetSeed(rctx, (unsigned long) (tval.tv_usec | getpid()));
    PetscRandomSeed(rctx);

    VecSetRandom(*(Vec*) x->ptr, rctx);
    PetscRandomDestroy(rctx);
#else
    (void) x;
#endif
}

void petsc_mat_duplicate_proc(adj_matrix matin, adj_matrix *matout) 
{
    /* Duplicate a matrix */
#ifdef HAVE_PETSC
    Mat *new_mat=(Mat*) malloc(sizeof(Mat));
    MatDuplicate(*(Mat*) matin.ptr, MAT_DO_NOT_COPY_VALUES, new_mat);
    MatZeroEntries(*new_mat);
    matout->ptr = (adj_matrix *) new_mat;
#else
    (void) matin;
    (void) matout;
#endif
}

void petsc_solve_proc(adj_variable var, adj_matrix mat, adj_vector rhs, adj_vector *soln) 
{
    /*************************************************/
    /*  Solve mat*soln=rhs using a direct LU solver  */
    /* You might want to implement your own callback */
    /* with different solver options.                */
    /*************************************************/
#ifdef HAVE_PETSC
    KSP            ksp; /* linear solver context */ 
    PC             pc;  /* preconditioner context */
    PetscTruth assembled;
   
    MatAssembled(*(Mat*) &mat, &assembled);
    if (!assembled)
      MatAssemblyEnd(petsc_mat_from_adj_matrix(mat), MAT_FINAL_ASSEMBLY);

    /* Create the output vector */
    Vec *sol_vec=(Vec*) malloc(sizeof(Vec));
    MatGetVecs(petsc_mat_from_adj_matrix(mat), sol_vec, NULL);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, petsc_mat_from_adj_matrix(mat), petsc_mat_from_adj_matrix(mat), DIFFERENT_NONZERO_PATTERN);

    KSPGetPC(ksp, &pc);
    KSPSetType(ksp, KSPPREONLY);
    PCSetType(pc, PCLU);
    KSPSetTolerances(ksp, 1.e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
    
    KSPSolve(ksp,petsc_vec_from_adj_vector(rhs), *sol_vec);
    *soln=petsc_vec_to_adj_vector(sol_vec);
#else
    (void) mat;
    (void) rhs;
    (void) soln;
#endif
    (void) var;
}

void petsc_mat_axpy_proc(adj_matrix *Y, adj_scalar alpha, adj_matrix X)
{
    /* Computes Y = alpha*X + Y. */
#ifdef HAVE_PETSC
    MatAXPY(*(Mat*) Y->ptr, alpha, *(Mat*) X.ptr, SAME_NONZERO_PATTERN);
#else
    (void) Y;
    (void) alpha;
    (void) X;
#endif
}

void petsc_mat_destroy_proc(adj_matrix *mat)
{
    /* Frees space taken by a matrix. */
#ifdef HAVE_PETSC
    MatDestroy(*(Mat*) mat->ptr);
    free(mat->ptr);
#else
    (void) mat;
#endif
}

void petsc_vec_getsize_proc(adj_vector vec, int *sz)
{
#ifdef HAVE_PETSC
  VecGetSize(*(Vec*) vec.ptr, sz);
#else
  (void) vec;
  (void) sz;
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
#else
  (void) vec;
  (void) scalars;
#endif
}

void petsc_vec_divide_proc(adj_vector *numerator, adj_vector denominator)
{
#ifdef HAVE_PETSC
  VecPointwiseDivide(*(Vec*) numerator->ptr, *(Vec*) numerator->ptr, *(Vec*) denominator.ptr);
#else
  (void) numerator;
  (void) denominator;
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

adj_matrix petsc_mat_to_adj_matrix(Mat* v)
{
  adj_matrix vv;
  vv.ptr = (void*)v;
  vv.klass = 0;
  return vv;
}

Mat petsc_mat_from_adj_matrix(adj_matrix vv)
{
  Mat v;
  v = *(Mat*) vv.ptr;
  return v;
}
#endif
