#include "libadjoint/adj_gst.h"
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))

int adj_compute_gst(adj_adjointer* adjointer, adj_variable ic, adj_matrix* ic_norm, adj_variable final, adj_matrix* final_norm, int nrv, adj_gst* gst_handle, int* ncv, int which)
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
  return ADJ_ERR_SLEPC_ERROR;
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

  strncpy(adj_error_msg, "Need the ADJ_MAT_ACTION_CB data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (final_norm != NULL && adjointer->callbacks.mat_action == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  strncpy(adj_error_msg, "Need the ADJ_VEC_DOT_PRODUCT_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (final_norm != NULL && adjointer->callbacks.vec_dot_product == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  gst_data = (adj_gst_data*) malloc(sizeof(adj_gst_data));
  gst_data->adjointer = adjointer;
  gst_data->ic = ic;
  gst_data->ic_norm = ic_norm;
  gst_data->final = final;
  gst_data->final_norm = final_norm;
  gst_data->multiplications = 0;

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
  EPSSetProblemType(*eps, EPS_NHEP);
  EPSSetType(*eps, EPSKRYLOVSCHUR);

  nwv = 3*nrv;
  if (nwv > global_ic_dof)
    nwv = PETSC_DECIDE;

  EPSSetDimensions(*eps, nrv, PETSC_DECIDE, PETSC_DECIDE);
  EPSSetWhichEigenpairs(*eps, which);

  ierr = EPSSetFromOptions(*eps);

  /* ierr = EPSView(*eps, PETSC_VIEWER_STDOUT_WORLD); */
  ierr = EPSSolve(*eps);

  printf("GST calculation took %d multiplications of L^*L.\n", gst_data->multiplications);

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
  return ADJ_ERR_SLEPC_ERROR;
#else

  int ierr;
  EPS* eps;
  adj_scalar ssigma;
  adj_scalar ssigma_complex;
  Vec v_vec;
  adj_adjointer* adjointer;

  eps = (EPS*) gst_handle->eps_handle;
  adj_gst_data* gst_data = (adj_gst_data*) gst_handle->gst_data;
  adjointer = gst_data->adjointer;

  if (u != NULL || v != NULL || sigma != NULL) /* we need to pull the eigenfunction */
  {
    Mat A;

    /* Shut the compiler up about uninitialised variables */
    EPSGetOperators(*( (EPS*) gst_handle->eps_handle ), &A, PETSC_NULL);
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 5 && PETSC_VERSION_RELEASE == 1
    MatGetVecs(A, &v_vec, PETSC_NULL);
#else
    MatCreateVecs(A, &v_vec, PETSC_NULL);
#endif

    ierr = EPSGetEigenpair(*eps, i, &ssigma, &ssigma_complex, v_vec, PETSC_NULL);
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

  if (u != NULL) /* this is the output perturbation, which we need to compute, alas */
  {
    adj_vector final_val;
    Vec u_vec; /* the vector that contains the tlm output */
    adj_scalar* u_arr;

#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 5 && PETSC_VERSION_RELEASE == 1
    MatGetVecs(gst_data->tlm_mat, PETSC_NULL, &u_vec);
#else
    MatCreateVecs(gst_data->tlm_mat, PETSC_NULL, &u_vec);
#endif
    MatMult(gst_data->tlm_mat, v_vec, u_vec); /* do the TLM solve */

    ierr = adj_get_variable_value(adjointer, gst_data->final, &final_val);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adjointer->callbacks.vec_duplicate(final_val, u);

    if (gst_data->final_norm == NULL)
    {
      VecNormalize(u_vec, PETSC_NULL);
    }
    else
    {
      /* Get its norm wrt final_norm, then call VecScale. */
      /* we'll use u as a temporary storage variable here. */
      adj_vector Xu;
      adj_scalar inner;
      adj_scalar norm;

      ierr = VecGetArray(u_vec, &u_arr);
      adjointer->callbacks.vec_set_values(u, u_arr);
      ierr = VecRestoreArray(u_vec, &u_arr);

      /* Now u is an adj_vector containing the unnormalized values */
      adjointer->callbacks.vec_duplicate(*u, &Xu);
      adjointer->callbacks.mat_action(*gst_data->final_norm, *u, &Xu);

      /* Now Xu contains the product X.u. We want to inner that with u
         to get the actual norm */
      adjointer->callbacks.vec_dot_product(*u, Xu, &inner);
      adjointer->callbacks.vec_destroy(&Xu);

      /* Now finally scale u_vec by the norm */
      norm = sqrt(inner);
      VecScale(u_vec, 1.0/norm);
    }

    ierr = VecGetArray(u_vec, &u_arr);
    adjointer->callbacks.vec_set_values(u, u_arr);
    ierr = VecRestoreArray(u_vec, &u_arr);

    ierr = VecDestroy(&u_vec);
  }

  if (v != NULL) /* this is the input perturbation, which we already have, handily */
  {
    adj_vector ic_val;
    adj_scalar* v_arr;

    ierr = adj_get_variable_value(adjointer, gst_data->ic, &ic_val);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adjointer->callbacks.vec_duplicate(ic_val, v);

    ierr = VecGetArray(v_vec, &v_arr);
    adjointer->callbacks.vec_set_values(v, v_arr);
    ierr = VecRestoreArray(v_vec, &v_arr);

    if (gst_data->ic_norm != NULL) /* need to normalise */
    {
      adj_vector Xv;
      adj_scalar inner;
      adj_scalar norm;

      adjointer->callbacks.vec_duplicate(*v, &Xv);
      adjointer->callbacks.mat_action(*gst_data->ic_norm, *v, &Xv);
      adjointer->callbacks.vec_dot_product(*v, Xv, &inner);
      adjointer->callbacks.vec_destroy(&Xv);

      /* Now finally scale v_vec by the norm */
      norm = sqrt(inner);
      VecScale(v_vec, 1.0/norm);

      ierr = VecGetArray(v_vec, &v_arr);
      adjointer->callbacks.vec_set_values(v, v_arr);
      ierr = VecRestoreArray(v_vec, &v_arr);
    }

  }

  if (error != NULL)
  {
    ierr = EPSComputeError(*eps, i, EPS_ERROR_RELATIVE, error);
    if (ierr != 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSComputeRelativeError (ierr == %d).", ierr);
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
  }

  if (sigma != NULL)
  {
    if (abs(ssigma) < DBL_EPSILON && ssigma < 0.0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc returned a negative number as a growth rate (error in SLEPc)?");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
    else
    {
      *sigma = sqrt(ssigma);
    }

    if (isnan(*sigma))
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc returned NaN as a growth rate.");
      return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
    }
  }

  ierr = VecDestroy(&v_vec);

  return ADJ_OK;

#endif
}

int adj_destroy_gst(adj_gst* gst_handle)
{
#ifndef HAVE_SLEPC
  (void) gst_handle;

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to destroy EPS objects, you need to compile with SLEPc support.");
  gst_handle = (adj_gst*) NULL;
  return ADJ_ERR_SLEPC_ERROR;
#else
  int ierr;
  adj_gst_data* gst_data;
#if PETSC_VERSION_MINOR > 1
  ierr = EPSDestroy(( (EPS*) gst_handle->eps_handle ));
#else
  ierr = EPSDestroy(*( (EPS*) gst_handle->eps_handle ));
#endif
  if (ierr != 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "SLEPc error from EPSDestroy (ierr == %d).", ierr);
    return adj_chkierr_auto(ADJ_ERR_SLEPC_ERROR);
  }
  free((EPS*) gst_handle->eps_handle);
  gst_data = (adj_gst_data*) gst_handle->gst_data;
#if PETSC_VERSION_MINOR > 1
  MatDestroy(&gst_data->tlm_mat);
#else
  MatDestroy(gst_data->tlm_mat);
#endif
  free((adj_gst_data*) gst_handle->gst_data);

  return ADJ_OK;
#endif
}

#ifdef HAVE_SLEPC

#undef __FUNCT__
#define __FUNCT__ "tlm_solve"
PetscErrorCode tlm_solve(Mat A, Vec x, Vec y)
{
  adj_gst_data* gst_data;
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

  // Benchmark variables
  time_t start, end;
  int wall_time_used;
  static int total_time = 0;
  static int nb_tlm_solves = 0;

  PetscFunctionBegin;

  printf("Beginning tlm_solve\n");
  time(&start);

  ierr = MatShellGetContext(A, (void**) &gst_data); CHKERRQ(ierr);
  adjointer = gst_data->adjointer;
  ierr = adj_equation_count(adjointer, &equation_count);

  return_flag = ADJ_FALSE;

  for (equation = 0; equation < equation_count; equation++)
  {
    ierr = adj_get_tlm_equation(adjointer, equation, "GSTNullTLM", &lhs, &rhs, &tlm_var);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

    ierr = adj_create_variable(tlm_var.name, tlm_var.timestep, tlm_var.iteration, ADJ_FALSE, &fwd_var);
    if (adj_variable_equal(&gst_data->ic, &fwd_var, 1))
    {
      /* fetch the vector from our input PETSc Vec, stuff it into rhs_tmp */
      adjointer->callbacks.vec_duplicate(rhs, &rhs_tmp);

      ierr = VecGetArrayRead(x, (const PetscScalar**) &px); CHKERRQ(ierr);
      adjointer->callbacks.vec_set_values(&rhs_tmp, px);
      ierr = VecRestoreArrayRead(x, (const PetscScalar**) &px); CHKERRQ(ierr);

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

    if (adj_variable_equal(&gst_data->final, &fwd_var, 1))
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
      time(&end);
      wall_time_used =  end - start;
      total_time += wall_time_used;
      nb_tlm_solves += 1;
      printf("Ending tlm_solve. Run time: %i s. Average run time: %f s (averaged over %i tlm_solve's)\n", wall_time_used, ((double) total_time)/nb_tlm_solves, nb_tlm_solves);
      PetscFunctionReturn(0);
    }
  }

  time(&end);
  wall_time_used = end - start;
  total_time += wall_time_used;
  nb_tlm_solves += 1;
  printf("Ending tlm_solve. Run time: %i s. Average run time: %f s (averaged over %i tlm_solve's)\n", wall_time_used, ((double) total_time)/nb_tlm_solves, nb_tlm_solves);
  PetscFunctionReturn(1);
}

#undef __FUNCT__
#define __FUNCT__ "adj_solve"
PetscErrorCode adj_solve(Mat A, Vec x, Vec y)
{
  adj_gst_data* gst_data;
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

  // Benchmark variables
  time_t start, end;
  int wall_time_used;
  static int total_time = 0;
  static int nb_tlm_solves = 0;

  PetscFunctionBegin;

  printf("Beginning adj_solve\n");
  time(&start);

  ierr = MatShellGetContext(A, (void**) &gst_data); CHKERRQ(ierr);
  adjointer = gst_data->adjointer;
  ierr = adj_equation_count(adjointer, &equation_count);

  return_flag = ADJ_FALSE;

  for (equation = equation_count - 1; equation >= 0; equation--)
  {
    ierr = adj_get_adjoint_equation(adjointer, equation, "GSTNullADM", &lhs, &rhs, &adj_var);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

    ierr = adj_create_variable(adj_var.name, adj_var.timestep, adj_var.iteration, ADJ_FALSE, &fwd_var);
    if (adj_variable_equal(&gst_data->final, &fwd_var, 1))
    {
      /* fetch the vector from our input PETSc Vec, stuff it into rhs_tmp */
      adjointer->callbacks.vec_duplicate(rhs, &rhs_tmp);

      ierr = VecGetArrayRead(x, (const PetscScalar**) &px); CHKERRQ(ierr);
      adjointer->callbacks.vec_set_values(&rhs_tmp, px);
      ierr = VecRestoreArrayRead(x, (const PetscScalar**) &px); CHKERRQ(ierr);

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

    if (adj_variable_equal(&gst_data->ic, &fwd_var, 1))
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
      time(&end);
      wall_time_used =  end - start;
      total_time += wall_time_used;
      nb_tlm_solves += 1;
      printf("Ending adj_solve. Run time: %i s. Average run time: %f s (averaged over %i adj_solve's)\n", wall_time_used, ((double) total_time)/nb_tlm_solves, nb_tlm_solves);
      PetscFunctionReturn(0);
    }
  }

  time(&end);
  wall_time_used =  end - start;
  total_time += wall_time_used;
  nb_tlm_solves += 1;
  printf("Ending adj_solve. Run time: %i s. Average run time: %f s (averaged over %i adj_solve's)\n", wall_time_used, ((double) total_time)/nb_tlm_solves, nb_tlm_solves);
  PetscFunctionReturn(1);
}

#undef __FUNCT__
#define __FUNCT__ "gst_mult"
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
  gst_data->multiplications++;

  /* Multiply by L .. */
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 5 && PETSC_VERSION_RELEASE == 1
  ierr = MatGetVecs(tlm_mat, PETSC_NULL, &Lx);   CHKERRQ(ierr);
#else
  ierr = MatCreateVecs(tlm_mat, PETSC_NULL, &Lx);   CHKERRQ(ierr);
#endif
  ierr = MatMult(tlm_mat, x, Lx);                CHKERRQ(ierr);

  /* Then take the final norm */
  if (gst_data->final_norm != NULL)
  {
    adj_scalar* Lx_array; /* to do some annoying shuffling from PETSc Vecs -> adj_vectors */
    adj_scalar* XLx_array;
    adj_vector Lx_vector;
    adj_vector XLx_vector;
    adj_vector final_val;

    /* Do the necessary allocations */
    ierr = VecDuplicate(Lx, &XLx);                 CHKERRQ(ierr);
    ierr = VecGetArray(Lx, &Lx_array);             CHKERRQ(ierr);
    ierr = VecGetArray(XLx, &XLx_array);           CHKERRQ(ierr);

    ierr = adj_get_variable_value(adjointer, gst_data->final, &final_val);
    adjointer->callbacks.vec_duplicate(final_val, &Lx_vector);
    adjointer->callbacks.vec_duplicate(final_val, &XLx_vector);

    /* OK. Stuff the values from Lx into Lx_vector. */
    adjointer->callbacks.vec_set_values(&Lx_vector, Lx_array);
    /* Now compute the action of the matrix. (Sets XLx_vector)*/
    adjointer->callbacks.mat_action(*gst_data->final_norm, Lx_vector, &XLx_vector);
    /* Now fetch the values from XLx_vector into XLx_array. */
    adjointer->callbacks.vec_get_values(XLx_vector, &XLx_array);

    /* Now clean up */
    ierr = VecRestoreArray(Lx, &Lx_array);         CHKERRQ(ierr);
    ierr = VecRestoreArray(XLx, &XLx_array);       CHKERRQ(ierr);

    adjointer->callbacks.vec_destroy(&Lx_vector);
    adjointer->callbacks.vec_destroy(&XLx_vector);

    ierr = VecDestroy(&Lx);                         CHKERRQ(ierr);

    /* Now XLx contains the action of the norm matrix, and everything
       else is destroyed. */
  }
  else
  {
    ierr = VecDuplicate(Lx, &XLx);                 CHKERRQ(ierr);
    ierr = VecCopy(Lx, XLx);                       CHKERRQ(ierr);
    ierr = VecDestroy(&Lx);                         CHKERRQ(ierr);
  }

  /* Now multiply by L^* .. */
#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 5 && PETSC_VERSION_RELEASE == 1
  ierr = MatGetVecs(tlm_mat, &LXLx, PETSC_NULL); CHKERRQ(ierr);
#else
  ierr = MatCreateVecs(tlm_mat, &LXLx, PETSC_NULL); CHKERRQ(ierr);
#endif
  ierr = MatMultTranspose(tlm_mat, XLx, LXLx);   CHKERRQ(ierr);
  ierr = VecDestroy(&XLx);

  /* Now take the initial norm */
  if (gst_data->ic_norm == NULL)
  {
    ierr = VecCopy(LXLx, y);                       CHKERRQ(ierr);
  }
  else
  {
    /* We need to do y = ic_norm^{-1} . LXLx */
    adj_vector y_vec;
    adj_vector LXLx_vec;
    adj_vector ic_val;
    adj_scalar* LXLx_array;
    adj_scalar* y_array;

    ierr = adj_get_variable_value(adjointer, gst_data->ic, &ic_val);
    adjointer->callbacks.vec_duplicate(ic_val, &y_vec);
    adjointer->callbacks.vec_duplicate(ic_val, &LXLx_vec);

    /* Set LXLx_vec from the PETSc array */
    ierr = VecGetArray(LXLx, &LXLx_array);         CHKERRQ(ierr);
    adjointer->callbacks.vec_set_values(&LXLx_vec, LXLx_array);
    ierr = VecRestoreArray(LXLx, &LXLx_array);     CHKERRQ(ierr);

    /* Now do the solve */
    adjointer->callbacks.solve(gst_data->ic, *gst_data->ic_norm, LXLx_vec, &y_vec);

    /* Now set the values of y */
    ierr = VecGetArray(y, &y_array);               CHKERRQ(ierr);
    adjointer->callbacks.vec_get_values(y_vec, &y_array);
    ierr = VecRestoreArray(y, &y_array);           CHKERRQ(ierr);
  }
  ierr = VecDestroy(&LXLx);

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
