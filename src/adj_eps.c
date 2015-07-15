#include "libadjoint/adj_eps.h"
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

#ifdef HAVE_SLEPC
#undef __FUNCT__
#define __FUNCT__ "eps_mult"
PetscErrorCode eps_mult(Mat A, Vec x, Vec y)
{
  adj_eps_data* eps_data;
  adj_adjointer* adjointer;
  adj_matrix matrix;
  adj_scalar* input_arr;
  adj_scalar* output_arr;
  adj_vector work_input;
  adj_vector work_output;
  int ierr;


  PetscFunctionBegin;

  ierr = MatShellGetContext(A, (void**) &eps_data); CHKERRQ(ierr);
  adjointer = eps_data->adjointer;
  matrix = eps_data->matrix;
  eps_data->multiplications++;
  printf("Beginning matrix action %d.\n", eps_data->multiplications-1);

  adjointer->callbacks.vec_duplicate(eps_data->input, &work_input);
  ierr = VecGetArrayRead(x, &input_arr); CHKERRQ(ierr);
  adjointer->callbacks.vec_set_values(&work_input, input_arr);
  ierr = VecRestoreArrayRead(x, &input_arr); CHKERRQ(ierr);

  adjointer->callbacks.vec_duplicate(eps_data->output, &work_output);

  adjointer->callbacks.mat_action(matrix, work_input, &work_output);

  ierr = VecGetArray(y, &output_arr); CHKERRQ(ierr);
  adjointer->callbacks.vec_get_values(work_output, &output_arr);
  ierr = VecRestoreArray(y, &output_arr); CHKERRQ(ierr);

  adjointer->callbacks.vec_destroy(&work_input);
  adjointer->callbacks.vec_destroy(&work_output);
  printf("Matrix action %d completed.\n", eps_data->multiplications-1);

  PetscFunctionReturn(0);
}

#endif /* HAVE_SLEPC */

int adj_compute_eps(adj_adjointer* adjointer, adj_matrix matrix, adj_eps_options options, adj_eps* eps_handle, int* nconverged)
{
#ifndef HAVE_SLEPC
  (void) adjointer;
  (void) eps_handle;
  (void) options;
  (void) matrix;
  (void) nconverged;

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to perform an eigendecomposition, you need to compile with SLEPc support.");
  return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
#else
  EPS *eps;
  Mat eps_mat;
  adj_eps_data* eps_data;

  int ierr;
  int input_dof, output_dof;
  int global_input_dof, global_output_dof;
  int nwv; /* number of work vectors */

  /* Check for the required callbacks */
  strncpy(adj_error_msg, "Need the ADJ_SOLVE_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.solve == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  strncpy(adj_error_msg, "Need the ADJ_VEC_AXPY_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.vec_axpy == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  strncpy(adj_error_msg, "Need the ADJ_VEC_DUPLICATE_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.vec_duplicate == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  strncpy(adj_error_msg, "Need the ADJ_VEC_DESTROY_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.vec_destroy == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  strncpy(adj_error_msg, "Need the ADJ_MAT_DESTROY_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.mat_destroy == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  strncpy(adj_error_msg, "Need the ADJ_VEC_GET_SIZE_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.vec_get_size == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  strncpy(adj_error_msg, "Need the ADJ_VEC_GET_VALUES_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.vec_get_values == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  strncpy(adj_error_msg, "Need the ADJ_VEC_SET_VALUES_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.vec_set_values == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  strncpy(adj_error_msg, "Need the ADJ_MAT_ACTION_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.mat_action == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  eps_data = (adj_eps_data*) malloc(sizeof(adj_eps_data));
  eps_data->adjointer = adjointer;
  eps_data->matrix = matrix;
  eps_data->input = options.input;
  eps_data->output = options.output;
  eps_data->multiplications = 0;

  adjointer->callbacks.vec_get_size(eps_data->input, &input_dof);
  adjointer->callbacks.vec_get_size(eps_data->output, &output_dof);

  ierr = MatCreateShell(PETSC_COMM_WORLD, output_dof, input_dof, PETSC_DETERMINE, PETSC_DETERMINE, (void*) eps_data, &eps_mat);
  ierr = MatShellSetOperation(eps_mat, MATOP_MULT, (void(*)(void)) eps_mult);

  MatGetSize(eps_mat, &global_output_dof, &global_input_dof);

  if (options.neigenpairs > min(global_input_dof, global_output_dof))
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Cannot request %d vectors when the matrix is %d x %d.", options.neigenpairs, global_output_dof, global_input_dof);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  SlepcInitialize(0, PETSC_NULL, PETSC_NULL, PETSC_NULL);

  eps = (EPS*) malloc(sizeof(EPS));
  EPSCreate(PETSC_COMM_WORLD, eps);
  EPSSetOperators(*eps, eps_mat, PETSC_NULL);

  EPSSetProblemType(*eps, (EPSProblemType) options.type);
  EPSSetType(*eps, options.method);

  nwv = 3*options.neigenpairs;
  if (nwv > global_input_dof) nwv = PETSC_DECIDE;

  EPSSetDimensions(*eps, options.neigenpairs, PETSC_DECIDE, PETSC_DECIDE);
#if SLEPC_VERSION_MAJOR > 3 || (SLEPC_VERSION_MAJOR == 3 && SLEPC_VERSION_MINOR >= 1)
  if (options.monitor) EPSMonitorSet(*eps, EPSMonitorAll, PETSC_NULL, PETSC_NULL);
#endif
  EPSSetWhichEigenpairs(*eps, (EPSWhich) options.which);

  if (options.monitor) EPSView(*eps, PETSC_VIEWER_STDOUT_WORLD);
  ierr = EPSSolve(*eps);

  if (options.monitor) printf("Eigenvalue calculation took %d multiplications.\n", eps_data->multiplications);

  eps_handle->eps_handle = eps;
  eps_handle->eps_data = eps_data;

  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSSolve (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  ierr = EPSGetConverged(*eps, nconverged);
  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSGetConverged (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  return ADJ_OK;

#endif
}

int adj_get_eps(adj_eps* eps_handle, int i, adj_scalar* sigma_re, adj_scalar* sigma_im, adj_vector* u_re, adj_vector* u_im)
{
#ifndef HAVE_SLEPC
  (void) i;
  (void) sigma_re;
  (void) sigma_im;
  (void) u_re;
  (void) u_im;
  (void) eps_handle;

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to fetch the eigenvalue analysis, you need to compile with SLEPc support.");
  eps_handle = (adj_eps*) NULL;
  return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
#else

  int ierr;
  EPS* eps;
  adj_scalar ssigma_re;
  adj_scalar ssigma_im;
  Vec u_vec_re;
  Vec u_vec_im;
  adj_adjointer* adjointer;
  Mat A;
  adj_eps_data* eps_data;
  adj_scalar norm;

  eps = (EPS*) eps_handle->eps_handle;
  eps_data = (adj_eps_data*) eps_handle->eps_data;
  adjointer = eps_data->adjointer;

  if (sigma_re == NULL && u_re == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Must ask for at least one of the eigenvalue or eigenvector.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  /* Shut the compiler up about uninitialised variables */
  EPSGetOperators(*( (EPS*) eps_handle->eps_handle ), &A, PETSC_NULL);
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 5 && PETSC_VERSION_RELEASE == 1
  MatGetVecs(A, &u_vec_re, PETSC_NULL);
#else
  MatCreateVecs(A, &u_vec_re, PETSC_NULL);
#endif
  VecDuplicate(u_vec_re, &u_vec_im);

  ierr = EPSGetEigenpair(*eps, i, &ssigma_re, &ssigma_im, u_vec_re, u_vec_im);
  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSGetEigenpair (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  if (ssigma_im != 0 && sigma_im == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc returned a complex eigenvalue in EPSGetEigenpair, but you passed sigma_im == NULL.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  VecNorm(u_vec_im, NORM_INFINITY, &norm);
  if (norm != 0 && u_im == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc returned a complex eigenvector in EPSGetEigenpair, but you passed u_im == NULL.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }
  if (isnan(ssigma_re) || isnan(ssigma_im))
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc returned NaN as an eigenvalue.");
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  /* After all that checking, now set the user's requested outputs */
  if (sigma_re) *sigma_re = ssigma_re;
  if (sigma_im) *sigma_im = ssigma_im;
  if (u_re != NULL)
  {
    adj_scalar* u_arr;

    adjointer->callbacks.vec_duplicate(eps_data->input, u_re);

    ierr = VecGetArray(u_vec_re, &u_arr);
    adjointer->callbacks.vec_set_values(u_re, u_arr);
    ierr = VecRestoreArray(u_vec_re, &u_arr);
  }

  ierr = VecDestroy(&u_vec_re);

  if (u_im != NULL)
  {
    adj_scalar* u_arr;

    adjointer->callbacks.vec_duplicate(eps_data->input, u_im);

    ierr = VecGetArray(u_vec_im, &u_arr);
    adjointer->callbacks.vec_set_values(u_im, u_arr);
    ierr = VecRestoreArray(u_vec_im, &u_arr);
  }

  ierr = VecDestroy(&u_vec_im);

  return ADJ_OK;

#endif
}

int adj_destroy_eps(adj_eps* eps_handle)
{
#ifndef HAVE_SLEPC
  (void) eps_handle;

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to destroy EPS objects, you need to compile with SLEPc support.");
  eps_handle = (adj_eps*) NULL;
  return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
#else
  int ierr;
#if PETSC_VERSION_MINOR > 1
  ierr = EPSDestroy(( (EPS*) eps_handle->eps_handle ));
#else
  ierr = EPSDestroy(*( (EPS*) eps_handle->eps_handle ));
#endif

  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSDestroy (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  free((adj_eps_data*) eps_handle->eps_data);

  return ADJ_OK;
#endif
}

