#ifndef ADJ_PREDICTABILITY_H
#define ADJ_PREDICTABILITY_H

#include "adj_data_structures.h"
#include "adj_error_handling.h"
#include "adj_adjointer_routines.h"
#include "adj_core.h"

#include <math.h>
#include <float.h>

#ifdef HAVE_SLEPC
#include "slepceps.h"
#endif

typedef struct
{
  void* eps_handle;
  void* gst_data;
} adj_gst;

int adj_compute_gst(adj_adjointer* adjointer, adj_variable ic, adj_matrix* ic_norm, adj_variable final, adj_matrix* final_norm, int nrv, adj_gst* gst_handle, int* ncv);
int adj_get_gst(adj_gst* gst_handle, int i, adj_scalar* sigma, adj_vector* u, adj_vector* v, adj_scalar* error);
int adj_destroy_gst(adj_gst* gst_handle);

#ifndef ADJ_HIDE_FROM_USER
#ifdef HAVE_SLEPC
typedef struct
{
  adj_adjointer* adjointer;
  adj_variable ic;
  adj_matrix* ic_norm;
  adj_variable final;
  adj_matrix* final_norm;
  Mat gst_mat;
  Mat tlm_mat;
} adj_gst_data;

PetscErrorCode tlm_solve(Mat A, Vec x, Vec y);
PetscErrorCode adj_solve(Mat A, Vec x, Vec y);
PetscErrorCode gst_mult(Mat A, Vec x, Vec y);
void null_tlm_source(adj_adjointer* adjointer, int equation, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output, int* has_output);
void null_adj_source(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output);
#endif
#endif

#endif

