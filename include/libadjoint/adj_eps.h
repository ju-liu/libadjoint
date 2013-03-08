#ifndef ADJ_EPS_H
#define ADJ_EPS_H

#include "adj_data_structures.h"
#include "adj_error_handling.h"

typedef struct
{
  adj_vector input;    /* A model input vector; not modified */
  adj_vector output;   /* A model output vector; not modified */
  char* method;        /* What algorithm should be used? */
  int type;            /* What type of problem is it? */
  int which;           /* Which eigenvectors do you want? */
  int monitor;         /* Do you want to print out information about convergence? */
  int neigenpairs;     /* How many vectors do you want? */
} adj_eps_options;

#ifdef HAVE_SLEPC
#include "slepceps.h"
#endif /* HAVE_SLEPC */

typedef struct
{
  void* eps_handle;        /* A pointer to the EPS */
  void* eps_data;          /* Any data the eigendecomposition might need as context */
} adj_eps;

int adj_compute_eps(adj_adjointer* adjointer, adj_matrix matrix, adj_eps_options options, adj_eps* eps_handle, int* nconverged);
int adj_get_eps(adj_eps* eps_handle, int i, adj_scalar* sigma_re, adj_scalar* sigma_im, adj_vector* u_re, adj_vector* u_im);
int adj_destroy_eps(adj_eps* eps_handle);

#ifndef ADJ_HIDE_FROM_USER
typedef struct
{
  adj_adjointer* adjointer;
  adj_matrix matrix;
  adj_vector input;
  adj_vector output;
  int multiplications;
} adj_eps_data;

#ifdef HAVE_PETSC
PetscErrorCode eps_mult(Mat A, Vec x, Vec y);
#endif /* HAVE_PETSC */
#endif /* ADJ_HIDE_FROM_USER */

#endif
