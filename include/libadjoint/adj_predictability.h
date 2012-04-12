#ifndef ADJ_PREDICTABILITY_H
#define ADJ_PREDICTABILITY_H

#include "adj_data_structures.h"
#include "adj_error_handling.h"
#include "adj_adjointer_routines.h"
#include "adj_core.h"

#include <math.h>

#ifdef HAVE_SLEPC
#include "slepcsvd.h"
#endif

typedef struct
{
  void* svd_handle;
  void* svd_data;
} adj_svd;

int adj_compute_propagator_svd(adj_adjointer* adjointer, adj_variable ic, adj_variable final, int nsv, adj_svd* svd_handle, int* ncv);
int adj_get_svd(adj_svd* svd_handle, int i, adj_scalar* sigma, adj_vector* u, adj_vector* v, adj_scalar* error);
int adj_destroy_svd(adj_svd* svd_handle);

#ifndef ADJ_HIDE_FROM_USER
typedef struct
{
  adj_adjointer* adjointer;
  adj_variable ic;
  adj_variable final;
} adj_svd_data;

#ifdef HAVE_SLEPC
PetscErrorCode tlm_solve(Mat A, Vec x, Vec y);
PetscErrorCode adj_solve(Mat A, Vec x, Vec y);
void null_tlm_source(adj_adjointer* adjointer, int equation, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output, int* has_output);
void null_adj_source(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output);
#endif
#endif

#endif

