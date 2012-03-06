#include "libadjoint/adj_predictability.h"

int adj_compute_tlm_svd(adj_adjointer* adjointer, adj_variable ic, adj_variable final, int nsv, adj_svd* svd_handle, int* ncv)
{
#ifndef HAVE_SLEPC
  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to compute the SVD of your TLM, you need to compile with SLEPc support.");
  *svd_handle = NULL;
  return ADJ_ERR_INVALID_INPUTS;
#else
  SVD svd;
  Mat tlm_mat;
  adj_svd_data svd_data;
  int ierr;
  int N;

  svd_data.adjointer = adjointer;
  svd_data.ic = ic;
  svd_data.final = final;

  N = 10;

  ierr = MatCreateShell(PETSC_COMM_WORLD, N, N, PETSC_DETERMINE, PETSC_DETERMINE, (void*) &svd_data, &tlm_mat);
/*  ierr = MatShellSetOperation(tlm_mat, MATOP_MULT,(void(*)()) tlm_solve);
  ierr = MatShellSetOperation(tlm_mat, MATOP_MULT_TRANSPOSE, (void(*)()) adj_solve); */

  SVDCreate(PETSC_COMM_WORLD, &svd);
  SVDSetOperator(svd, tlm_mat);
  SVDSetTransposeMode(svd, SVD_TRANSPOSE_IMPLICIT);
  SVDSetType(svd, SVDCROSS);
  SVDSetDimensions(svd, nsv, PETSC_DECIDE, PETSC_DECIDE);

  ierr = SVDSolve(svd);

  svd_handle->svd_handle = &svd;
  SVDGetConverged(svd, ncv);

  if (ierr != 0)
    return ADJ_OK;
  else
  {
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

#endif
}

int adj_get_svd(adj_svd* svd_handle, int i, adj_scalar* sigma, adj_vector* u, adj_vector* v, adj_scalar* error)
{
#ifndef HAVE_SLEPC
  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to fetch SVDs, you need to compile with SLEPc support.");
  *svd_handle = NULL;
  return ADJ_ERR_INVALID_INPUTS;
#else

  int ierr;

  if (u != NULL && v != NULL)
  {
    Mat A;
    Vec u_vec;
    Vec v_vec;

    /* Shut the compiler up about uninitialised variables */
    SVDGetOperator((SVD) svd_handle->svd_handle, &A);
    MatGetVecs(A, &u_vec, &v_vec);

    ierr = SVDGetSingularTriplet((SVD) svd_handle->svd_handle, i, sigma, u_vec, v_vec);

    if (ierr != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDGetSingularTriplet.");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }

    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Fetching singular vectors from the SVD not implemented (yet).");
    return adj_chkierr_auto(ADJ_ERR_NOT_IMPLEMENTED);

  }
  else
  {
    ierr = SVDGetSingularTriplet((SVD) svd_handle->svd_handle, i, sigma, PETSC_NULL, PETSC_NULL);
    if (ierr != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDGetSingularTriplet.");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
  }

  if (error != NULL)
  {
    ierr = SVDComputeRelativeError((SVD) svd_handle->svd_handle, i, error);
    if (ierr != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDComputeRelativeError.");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
  }

  return ADJ_OK;

#endif
}

int adj_destroy_svd(adj_svd* svd_handle)
{
#ifndef HAVE_SLEPC
  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to destroy SVD objects, you need to compile with SLEPc support.");
  *svd_handle = NULL;
  return ADJ_ERR_INVALID_INPUTS;
#else
  int ierr;
  ierr = SVDDestroy((SVD) svd_handle->svd_handle);
  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDDestroy.");
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  return ADJ_OK;
#endif
}

#ifdef HAVE_SLEPC
PetscErrorCode tlm_solve(Mat A, Vec x, Vec y);
PetscErrorCode adj_solve(Mat A, Vec x, Vec y);
#endif
