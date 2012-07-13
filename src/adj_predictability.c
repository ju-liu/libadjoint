#include "libadjoint/adj_predictability.h"
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

int adj_compute_propagator_svd(adj_adjointer* adjointer, adj_variable ic, adj_variable final, int nsv, adj_svd* svd_handle, int* ncv)
{
#ifndef HAVE_SLEPC
  (void) adjointer;
  (void) ic;
  (void) final;
  (void) nsv; /* number of requested vectors */
  (void) svd_handle;
  (void) ncv; /* number of converged vectors */

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to compute an SVD, you need to compile with SLEPc support.");
  svd_handle = (adj_svd*) NULL;
  return ADJ_ERR_INVALID_INPUTS;
#else
  SVD *svd;
  Mat tlm_mat;
  adj_gst_data* svd_data;
  int ierr;
  adj_vector ic_val;
  adj_vector final_val;
  int ic_dof, final_dof;
  int global_ic_dof, global_final_dof;
  int nwv; /* number of work vectors */

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

  svd_data = (adj_gst_data*) malloc(sizeof(adj_gst_data));
  svd_data->adjointer = adjointer;
  svd_data->ic = ic;
  svd_data->final = final;

  /* Register the dummy parameter sources/functionals -- we're going to be in charge of the RHS terms here */
  adj_register_parameter_source_callback(adjointer, "GSTNullTLM", null_tlm_source);
  adj_register_functional_derivative_callback(adjointer, "GSTNullADM", null_adj_source);

  ierr = adj_get_variable_value(adjointer, ic, &ic_val);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
  adjointer->callbacks.vec_get_size(ic_val, &ic_dof);

  ierr = adj_get_variable_value(adjointer, final, &final_val);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
  adjointer->callbacks.vec_get_size(final_val, &final_dof);

  ierr = MatCreateShell(PETSC_COMM_WORLD, final_dof, ic_dof, PETSC_DETERMINE, PETSC_DETERMINE, (void*) svd_data, &tlm_mat);
  ierr = MatShellSetOperation(tlm_mat, MATOP_MULT, (void(*)(void)) tlm_solve);
  ierr = MatShellSetOperation(tlm_mat, MATOP_MULT_TRANSPOSE, (void(*)(void)) adj_solve);

  MatGetSize(tlm_mat, &global_final_dof, &global_ic_dof);

  if (nsv > min(global_ic_dof, global_final_dof))
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Cannot request %d vectors when the matrix is %d x %d.", nsv, global_final_dof, global_ic_dof);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  SlepcInitialize(0, PETSC_NULL, PETSC_NULL, PETSC_NULL);

  svd = (SVD*) malloc(sizeof(SVD));
  SVDCreate(PETSC_COMM_WORLD, svd);
  SVDSetOperator(*svd, tlm_mat);
  SVDSetTransposeMode(*svd, SVD_TRANSPOSE_IMPLICIT);
  SVDSetType(*svd, SVDTRLANCZOS);
  SVDTRLanczosSetOneSide(*svd, PETSC_FALSE);

  nwv = 3*nsv;
  if (nwv > min(global_ic_dof, global_final_dof))
    nwv = PETSC_DECIDE;

  SVDSetDimensions(*svd, nsv, nwv, PETSC_DECIDE);
#if SLEPC_VERSION_MAJOR > 3 || (SLEPC_VERSION_MAJOR == 3 && SLEPC_VERSION_MINOR >= 1)
  SVDMonitorSet(*svd, SVDMonitorAll, PETSC_NULL, PETSC_NULL);
#endif

  ierr = SVDSolve(*svd);

  svd_handle->svd_handle = svd;
  svd_handle->svd_data = svd_data;

  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDSolve (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  ierr = SVDGetConverged(*svd, ncv);
  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDGetConverged (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  return ADJ_OK;

#endif
}

int adj_get_svd(adj_svd* svd_handle, int i, adj_scalar* sigma, adj_vector* u, adj_vector* v, adj_scalar* error)
{
#ifndef HAVE_SLEPC
  (void) i;
  (void) sigma;
  (void) u;
  (void) v;
  (void) error;
  (void) svd_handle;

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to fetch SVDs, you need to compile with SLEPc support.");
  svd_handle = (adj_svd*) NULL;
  return ADJ_ERR_INVALID_INPUTS;
#else

  int ierr;

  if (u != NULL && v != NULL)
  {
    Mat A;
    Vec u_vec;
    Vec v_vec;
    adj_vector ic_val;
    adj_vector final_val;
    adj_scalar* u_arr;
    adj_scalar* v_arr;
    adj_adjointer* adjointer;

    /* Shut the compiler up about uninitialised variables */
    SVDGetOperator(*( (SVD*) svd_handle->svd_handle ), &A);
    MatGetVecs(A, &v_vec, &u_vec);

    ierr = SVDGetSingularTriplet(*( (SVD*) svd_handle->svd_handle ), i, sigma, u_vec, v_vec);

    if (ierr != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDGetSingularTriplet (ierr == %d).", ierr);
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
    if (isnan(*sigma))
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc returned NaN as a singular value.");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }

    adjointer = ((adj_gst_data*) svd_handle->svd_data)->adjointer;

    ierr = adj_get_variable_value(adjointer, ((adj_gst_data*) svd_handle->svd_data)->ic, &ic_val);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adjointer->callbacks.vec_duplicate(ic_val, v);

    ierr = VecGetArray(v_vec, &v_arr);
    adjointer->callbacks.vec_set_values(v, v_arr);
    ierr = VecRestoreArray(v_vec, &v_arr);

    ierr = adj_get_variable_value(adjointer, ((adj_gst_data*) svd_handle->svd_data)->final, &final_val);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adjointer->callbacks.vec_duplicate(final_val, u);

    ierr = VecGetArray(u_vec, &u_arr);
    adjointer->callbacks.vec_set_values(u, u_arr);
    ierr = VecRestoreArray(u_vec, &u_arr);
  }
  else
  {
    ierr = SVDGetSingularTriplet(*( (SVD*) svd_handle->svd_handle ), i, sigma, PETSC_NULL, PETSC_NULL);
    if (ierr != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDGetSingularTriplet.");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
  }

  if (error != NULL)
  {
    ierr = SVDComputeRelativeError(*( (SVD*) svd_handle->svd_handle ), i, error);
    if (ierr != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDComputeRelativeError.");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
  }

  return ADJ_OK;

#endif
}

int adj_get_gst(adj_gst* gst_handle, int i, adj_scalar* sigma, adj_vector* u, adj_vector* v, adj_scalar* error)
{
#ifndef HAVE_SLEPC
  (void) i;
  (void) sigma;
  (void) u;
  (void) v;
  (void) error;
  (void) gst_handle;

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to fetch the GST analysis, you need to compile with SLEPc support.");
  gst_handle = (adj_gst*) NULL;
  return ADJ_ERR_INVALID_INPUTS;
#else

  int ierr;
  EPS* eps;
  adj_scalar ssigma;
  adj_scalar ssigma_complex;
  Vec v_vec;
  adj_adjointer* adjointer;
  
  eps = (EPS*) gst_handle->eps_handle;
  adjointer = ((adj_gst_data*) gst_handle->gst_data)->adjointer;

  if (u != NULL || v != NULL || sigma != NULL) /* we need to pull the eigenfunction */
  {
    Vec dummy;
    Mat A;

    /* Shut the compiler up about uninitialised variables */
    EPSGetOperators(*( (EPS*) gst_handle->eps_handle ), &A, PETSC_NULL);
    MatGetVecs(A, &v_vec, &dummy);

    ierr = EPSGetEigenpair(*eps, i, &ssigma, &ssigma_complex, v_vec, dummy);
    if (ierr != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSGetEigenpair (ierr == %d).", ierr);
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }

    if (ssigma_complex != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc returned a complex eigenvalue in EPSGetEigenpair.");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }

    if (isnan(ssigma))
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc returned NaN as an eigenvalue.");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
  }

  if (sigma != NULL)
  {
    *sigma = ssigma;
  }

  if (v != NULL) /* this is the input perturbation, which we already have, handily */
  {
    adj_vector ic_val;
    adj_scalar* v_arr;

    ierr = adj_get_variable_value(adjointer, ((adj_gst_data*) gst_handle->gst_data)->ic, &ic_val);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adjointer->callbacks.vec_duplicate(ic_val, v);

    ierr = VecGetArray(v_vec, &v_arr);
    adjointer->callbacks.vec_set_values(v, v_arr);
    ierr = VecRestoreArray(v_vec, &v_arr);
  }

  if (u != NULL) /* this is the output perturbation, which we need to compute, alas */
  {
    adj_vector final_val;
    Vec u_vec; /* the vector that contains the tlm output */
    adj_scalar* u_arr;

    MatGetVecs(((adj_gst_data*) gst_handle->gst_data)->tlm_mat, PETSC_NULL, &u_vec);
    MatMult(((adj_gst_data*) gst_handle->gst_data)->tlm_mat, v_vec, u_vec); /* do the TLM solve */

    ierr = adj_get_variable_value(adjointer, ((adj_gst_data*) gst_handle->gst_data)->final, &final_val);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adjointer->callbacks.vec_duplicate(final_val, u);

    ierr = VecGetArray(u_vec, &u_arr);
    adjointer->callbacks.vec_set_values(u, u_arr);
    ierr = VecRestoreArray(u_vec, &u_arr);
  }

  if (error != NULL)
  {
    ierr = EPSComputeRelativeError(*eps, i, error);
    if (ierr != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSComputeRelativeError (ierr == %d).", ierr);
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
  }

  return ADJ_OK;

#endif
}

int adj_destroy_svd(adj_svd* svd_handle)
{
#ifndef HAVE_SLEPC
  (void) svd_handle;

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to destroy SVD objects, you need to compile with SLEPc support.");
  svd_handle = (adj_svd*) NULL;
  return ADJ_ERR_INVALID_INPUTS;
#else
  int ierr;
  ierr = SVDDestroy(*( (SVD*) svd_handle->svd_handle ));
  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDDestroy.");
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }
  free((SVD*) svd_handle->svd_handle);
  free((adj_gst_data*) svd_handle->svd_data);

  return ADJ_OK;
#endif
}

int adj_compute_gst(adj_adjointer* adjointer, adj_variable ic, adj_matrix* ic_norm, adj_variable final, adj_matrix* final_norm, int nrv, adj_gst* gst_handle, int* ncv)
{
#ifndef HAVE_SLEPC
  (void) adjointer;
  (void) ic;
  (void) ic_norm;
  (void) final;
  (void) final_norm;
  (void) nrv; /* number of requested vectors */
  (void) gst_handle;
  (void) ncv; /* number of converged vectors */

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to perform generalised stability analysis, you need to compile with SLEPc support.");
  gst_handle = (adj_gst*) NULL;
  return ADJ_ERR_INVALID_INPUTS;
#else
  EPS *eps;
  Mat gst_mat;
  Mat tlm_mat;
  adj_gst_data* gst_data;

  int ierr;
  adj_vector ic_val;
  adj_vector final_val;
  int ic_dof, final_dof;
  int global_ic_dof, global_final_dof;
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

  gst_data = (adj_gst_data*) malloc(sizeof(adj_gst_data));
  gst_data->adjointer = adjointer;
  gst_data->ic = ic;
  gst_data->ic_norm = ic_norm;
  gst_data->final = final;
  gst_data->final_norm = final_norm;

  /* Register the dummy parameter sources/functionals -- we're going to be in charge of the RHS terms here */
  adj_register_parameter_source_callback(adjointer, "GSTNullTLM", null_tlm_source);
  adj_register_functional_derivative_callback(adjointer, "GSTNullADM", null_adj_source);

  ierr = adj_get_variable_value(adjointer, ic, &ic_val);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
  adjointer->callbacks.vec_get_size(ic_val, &ic_dof);

  ierr = adj_get_variable_value(adjointer, final, &final_val);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
  adjointer->callbacks.vec_get_size(final_val, &final_dof);

  ierr = MatCreateShell(PETSC_COMM_WORLD, ic_dof, ic_dof, PETSC_DETERMINE, PETSC_DETERMINE, (void*) gst_data, &gst_mat);
  ierr = MatShellSetOperation(gst_mat, MATOP_MULT, (void(*)(void)) gst_mult);

  ierr = MatCreateShell(PETSC_COMM_WORLD, final_dof, ic_dof, PETSC_DETERMINE, PETSC_DETERMINE, (void*) gst_data, &tlm_mat);
  ierr = MatShellSetOperation(tlm_mat, MATOP_MULT, (void(*)(void)) tlm_solve);
  ierr = MatShellSetOperation(tlm_mat, MATOP_MULT_TRANSPOSE, (void(*)(void)) adj_solve);
  gst_data->tlm_mat = tlm_mat;

  MatGetSize(tlm_mat, &global_final_dof, &global_ic_dof);

  if (nrv > min(global_ic_dof, global_final_dof))
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Cannot request %d vectors when the matrix is %d x %d.", nrv, global_final_dof, global_ic_dof);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  SlepcInitialize(0, PETSC_NULL, PETSC_NULL, PETSC_NULL);

  eps = (EPS*) malloc(sizeof(EPS));
  EPSCreate(PETSC_COMM_WORLD, eps);
  EPSSetOperators(*eps, gst_mat, PETSC_NULL);
  EPSSetProblemType(*eps, EPS_HEP);

  nwv = 3*nrv;
  if (nwv > global_ic_dof)
    nwv = PETSC_DECIDE;

  EPSSetDimensions(*eps, nrv, nwv, PETSC_DECIDE);
#if SLEPC_VERSION_MAJOR > 3 || (SLEPC_VERSION_MAJOR == 3 && SLEPC_VERSION_MINOR >= 1)
  EPSMonitorSet(*eps, EPSMonitorAll, PETSC_NULL, PETSC_NULL);
#endif

  ierr = EPSSolve(*eps);

  gst_handle->eps_handle = eps;
  gst_handle->gst_data = gst_data;

  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSSolve (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  ierr = EPSGetConverged(*eps, ncv);
  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSGetConverged (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  return ADJ_OK;

#endif
}

int adj_destroy_gst(adj_gst* gst_handle)
{
#ifndef HAVE_SLEPC
  (void) gst_handle;

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to destroy EPS objects, you need to compile with SLEPc support.");
  gst_handle = (adj_gst*) NULL;
  return ADJ_ERR_INVALID_INPUTS;
#else
  int ierr;
  adj_gst_data* gst_data;
  ierr = EPSDestroy(*( (EPS*) gst_handle->eps_handle ));
  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSDestroy (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }
  free((EPS*) gst_handle->eps_handle);
  gst_data = (adj_gst_data*) gst_handle->gst_data;
  MatDestroy(gst_data->tlm_mat);
  free((adj_gst_data*) gst_handle->gst_data);

  return ADJ_OK;
#endif
}

#ifdef HAVE_SLEPC
PetscErrorCode tlm_solve(Mat A, Vec x, Vec y)
{
  adj_gst_data* svd_data;
  adj_matrix lhs;
  adj_vector rhs;
  adj_vector rhs_tmp;
  adj_vector soln;
  adj_variable tlm_var;
  adj_variable fwd_var;
  adj_storage_data storage;
  int equation;
  int equation_count;
  int return_flag;
  adj_adjointer* adjointer;
  adj_scalar* px;
  adj_scalar* py;
  int ierr;

  PetscFunctionBegin;

  printf("Beginning tlm_solve\n");

  ierr = MatShellGetContext(A, (void**) &svd_data); CHKERRQ(ierr);
  adjointer = svd_data->adjointer;
  ierr = adj_equation_count(adjointer, &equation_count);

  return_flag = ADJ_FALSE;

  for (equation = 0; equation < equation_count; equation++)
  {
    ierr = adj_get_tlm_equation(adjointer, equation, "GSTNullTLM", &lhs, &rhs, &tlm_var);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

    ierr = adj_create_variable(tlm_var.name, tlm_var.timestep, tlm_var.iteration, ADJ_FALSE, &fwd_var);
    if (adj_variable_equal(&svd_data->ic, &fwd_var, 1))
    {
      /* fetch the vector from our input PETSc Vec, stuff it into rhs_tmp */
      adjointer->callbacks.vec_duplicate(rhs, &rhs_tmp);

      ierr = VecGetArray(x, &px); CHKERRQ(ierr);
      adjointer->callbacks.vec_set_values(&rhs_tmp, px);
      ierr = VecRestoreArray(x, &px); CHKERRQ(ierr);

      adjointer->callbacks.vec_axpy(&rhs, (adj_scalar) 1.0, rhs_tmp);
      adjointer->callbacks.vec_destroy(&rhs_tmp);
    }

    adjointer->callbacks.solve(tlm_var, lhs, rhs, &soln);
    adjointer->callbacks.vec_destroy(&rhs);
    adjointer->callbacks.mat_destroy(&lhs);

    ierr = adj_storage_memory_copy(soln, &storage);
    ierr = adj_storage_set_overwrite(&storage, ADJ_TRUE);
    ierr = adj_record_variable(adjointer, tlm_var, storage);

    ierr = adj_forget_tlm_values(adjointer, equation);

    if (adj_variable_equal(&svd_data->final, &fwd_var, 1))
    {
      /* fetch the value of soln, stuff it into our output PETSc Vec */
      ierr = VecGetArray(y, &py); CHKERRQ(ierr);
      adjointer->callbacks.vec_get_values(soln, &py);
      ierr = VecRestoreArray(y, &py); CHKERRQ(ierr);

      return_flag = ADJ_TRUE;
    }

    adjointer->callbacks.vec_destroy(&soln);

    if (return_flag) 
    {
      printf("Ending tlm_solve\n");
      PetscFunctionReturn(0);
    }
  }

  printf("Ending tlm_solve\n");
  PetscFunctionReturn(1);
}

PetscErrorCode adj_solve(Mat A, Vec x, Vec y)
{
  adj_gst_data* svd_data;
  adj_matrix lhs;
  adj_vector rhs;
  adj_vector rhs_tmp;
  adj_vector soln;
  adj_variable adj_var;
  adj_variable fwd_var;
  adj_storage_data storage;
  int equation;
  int equation_count;
  int return_flag;
  adj_adjointer* adjointer;
  adj_scalar* px;
  adj_scalar* py;
  int ierr;

  PetscFunctionBegin;

  printf("Beginning adj_solve\n");

  ierr = MatShellGetContext(A, (void**) &svd_data); CHKERRQ(ierr);
  adjointer = svd_data->adjointer;
  ierr = adj_equation_count(adjointer, &equation_count);

  return_flag = ADJ_FALSE;

  for (equation = equation_count - 1; equation >= 0; equation--)
  {
    ierr = adj_get_adjoint_equation(adjointer, equation, "GSTNullADM", &lhs, &rhs, &adj_var);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

    ierr = adj_create_variable(adj_var.name, adj_var.timestep, adj_var.iteration, ADJ_FALSE, &fwd_var);
    if (adj_variable_equal(&svd_data->final, &fwd_var, 1))
    {
      /* fetch the vector from our input PETSc Vec, stuff it into rhs_tmp */
      adjointer->callbacks.vec_duplicate(rhs, &rhs_tmp);

      ierr = VecGetArray(x, &px); CHKERRQ(ierr);
      adjointer->callbacks.vec_set_values(&rhs_tmp, px);
      ierr = VecRestoreArray(x, &px); CHKERRQ(ierr);

      adjointer->callbacks.vec_axpy(&rhs, (adj_scalar) 1.0, rhs_tmp);
      adjointer->callbacks.vec_destroy(&rhs_tmp);
    }

    adjointer->callbacks.solve(adj_var, lhs, rhs, &soln);
    adjointer->callbacks.vec_destroy(&rhs);
    adjointer->callbacks.mat_destroy(&lhs);

    ierr = adj_storage_memory_copy(soln, &storage);
    ierr = adj_storage_set_overwrite(&storage, ADJ_TRUE);
    ierr = adj_record_variable(adjointer, adj_var, storage);

    ierr = adj_forget_adjoint_values(adjointer, equation);

    if (adj_variable_equal(&svd_data->ic, &fwd_var, 1))
    {
      /* fetch the value of soln, stuff it into our output PETSc Vec */
      ierr = VecGetArray(y, &py); CHKERRQ(ierr);
      adjointer->callbacks.vec_get_values(soln, &py);
      ierr = VecRestoreArray(y, &py); CHKERRQ(ierr);

      return_flag = ADJ_TRUE;
    }

    adjointer->callbacks.vec_destroy(&soln);

    if (return_flag) 
    {
      printf("Ending adj_solve\n");
      PetscFunctionReturn(0);
    }
  }

  printf("Ending adj_solve\n");
  PetscFunctionReturn(1);
}

PetscErrorCode gst_mult(Mat A, Vec x, Vec y)
{
  adj_gst_data* gst_data;
  adj_adjointer* adjointer;
  
  Mat tlm_mat;
  Vec Lx;
  Vec XLx;
  Vec LXLx;

  int ierr;

  PetscFunctionBegin;

  ierr = MatShellGetContext(A, (void**) &gst_data); CHKERRQ(ierr);
  adjointer = gst_data->adjointer;
  tlm_mat = gst_data->tlm_mat;

  /* Multiply by L .. */
  ierr = MatGetVecs(tlm_mat, PETSC_NULL, &Lx);   CHKERRQ(ierr);
  ierr = MatMult(tlm_mat, x, Lx);                CHKERRQ(ierr);

  /* Then take the final norm */
  assert(gst_data->final_norm == NULL); /* for now */
  ierr = VecDuplicate(Lx, &XLx);                 CHKERRQ(ierr);
  ierr = VecCopy(Lx, XLx);                       CHKERRQ(ierr);
  ierr = VecDestroy(Lx);                         CHKERRQ(ierr);

  /* Now multiply by L^* .. */
  ierr = MatGetVecs(tlm_mat, &LXLx, PETSC_NULL); CHKERRQ(ierr);
  ierr = MatMultTranspose(tlm_mat, XLx, LXLx);   CHKERRQ(ierr);
  ierr = VecDestroy(XLx);

  /* Now take the initial norm */
  assert(gst_data->ic_norm == NULL); /* for now */
  ierr = VecCopy(LXLx, y);                       CHKERRQ(ierr);
  ierr = VecDestroy(LXLx);

  PetscFunctionReturn(0);
}

void null_tlm_source(adj_adjointer* adjointer, int equation, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output, int* has_output)
{
  (void) adjointer;
  (void) equation;
  (void) derivative;
  (void) ndepends;
  (void) variables;
  (void) dependencies;
  (void) name;
  (void) output;
  *has_output = ADJ_FALSE;
}

void null_adj_source(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output)
{
  (void) adjointer;
  (void) derivative;
  (void) ndepends;
  (void) variables;
  (void) dependencies;
  (void) name;
  (void) output;
}
#endif
