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
  adj_vector ic_val;
  adj_vector final_val;
  int ic_dof, final_dof;

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

  svd_data.adjointer = adjointer;
  svd_data.ic = ic;
  svd_data.final = final;

  /* Register the dummy parameter sources/functionals -- we're going to be in charge of the RHS terms here */
  adj_register_parameter_source_callback(adjointer, "SVDNullTLM", null_tlm_source);
  adj_register_functional_derivative_callback(adjointer, "SVDNullADM", null_adj_source);

  ierr = adj_get_variable_value(adjointer, ic, &ic_val);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
  adjointer->callbacks.vec_get_size(ic_val, &ic_dof);

  ierr = adj_get_variable_value(adjointer, final, &final_val);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
  adjointer->callbacks.vec_get_size(final_val, &final_dof);

  ierr = MatCreateShell(PETSC_COMM_WORLD, final_dof, ic_dof, PETSC_DETERMINE, PETSC_DETERMINE, (void*) &svd_data, &tlm_mat);
  ierr = MatShellSetOperation(tlm_mat, MATOP_MULT, (void(*)(void)) tlm_solve);
  ierr = MatShellSetOperation(tlm_mat, MATOP_MULT_TRANSPOSE, (void(*)(void)) adj_solve);

  SlepcInitialize(0, PETSC_NULL, PETSC_NULL, PETSC_NULL);

  SVDCreate(PETSC_COMM_WORLD, &svd);
  SVDSetOperator(svd, tlm_mat);
  SVDSetTransposeMode(svd, SVD_TRANSPOSE_IMPLICIT);
  SVDSetType(svd, SVDCROSS);
  SVDSetDimensions(svd, nsv, PETSC_DECIDE, PETSC_DECIDE);

  ierr = SVDSolve(svd);

  svd_handle->svd_handle = &svd;

  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDSolve.");
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  ierr = SVDGetConverged(svd, ncv);
  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from SVDGetConverged.");
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }

  return ADJ_OK;

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
PetscErrorCode tlm_solve(Mat A, Vec x, Vec y)
{
  adj_svd_data* svd_data;
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

  ierr = MatShellGetContext(A, (void**) &svd_data); CHKERRQ(ierr);
  adjointer = svd_data->adjointer;
  ierr = adj_equation_count(adjointer, &equation_count);

  return_flag = ADJ_FALSE;

  for (equation = 0; equation < equation_count; equation++)
  {
    ierr = adj_get_tlm_equation(adjointer, equation, "SVDNullTLM", &lhs, &rhs, &tlm_var);
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
      PetscFunctionReturn(0);
  }

  PetscFunctionReturn(1);
}

PetscErrorCode adj_solve(Mat A, Vec x, Vec y)
{
  adj_svd_data* svd_data;
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

  ierr = MatShellGetContext(A, (void**) &svd_data); CHKERRQ(ierr);
  adjointer = svd_data->adjointer;
  ierr = adj_equation_count(adjointer, &equation_count);

  return_flag = ADJ_FALSE;

  for (equation = equation_count - 1; equation >= 0; equation--)
  {
    ierr = adj_get_adjoint_equation(adjointer, equation, "SVDNullADM", &lhs, &rhs, &adj_var);
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
      PetscFunctionReturn(0);
  }

  PetscFunctionReturn(1);
}

void null_tlm_source(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output, int* has_output)
{
  (void) adjointer;
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
