#ifndef ADJ_PREDICTABILITY_H
#define ADJ_PREDICTABILITY_H

#include "adj_data_structures.h"
#include "adj_error_handling.h"

#ifdef HAVE_SLEPC
#include "slepcsvd.h"
#endif

typedef struct
{
  adj_adjointer* adjointer;
  adj_variable ic;
  adj_variable final;
} adj_svd_data;

typedef struct
{
  void* svd_handle;
} adj_svd;

int adj_compute_tlm_svd(adj_adjointer* adjointer, adj_variable ic, adj_variable final, int nsv, adj_svd* svd_handle, int* ncv);
int adj_get_svd(adj_svd* svd_handle, int i, adj_scalar* sigma, adj_vector* u, adj_vector* v, adj_scalar* error);
int adj_destroy_svd(adj_svd* svd_handle);

#ifdef HAVE_SLEPC
PetscErrorCode tlm_solve(Mat A, Vec x, Vec y);
PetscErrorCode adj_solve(Mat A, Vec x, Vec y);
#endif

#endif

