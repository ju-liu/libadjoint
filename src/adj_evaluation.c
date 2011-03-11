#include "libadjoint/adj_evaluation.h"

int adj_evaluate_block_action(adj_adjointer* adjointer, adj_block block, adj_vector input, adj_vector* output)
{
  int i, ierr;
  void (*block_action_func)(int, adj_variable*, adj_vector*, int, adj_scalar, adj_vector, void*, adj_vector*) = NULL;
  adj_vector* dependencies = NULL;
  int nb_variables = 0;
  adj_variable* variables = NULL;

  ierr = adj_find_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, block.name, (void (**)(void)) &block_action_func);
  if (ierr != ADJ_ERR_OK)
    return ierr;

  /* We have the right callback, so let's call it already */ 
  if (block.has_nonlinear_block) 
  { 
    /* we need to set up the dependencies */
    nb_variables = block.nonlinear_block.ndepends;
    variables = block.nonlinear_block.depends;
    dependencies = (adj_vector*) malloc(nb_variables * sizeof(adj_vector));

    for (i = 0; i < nb_variables; i++)
    {
      ierr = adj_get_variable_value(adjointer, variables[i], &dependencies[i]);
      if (ierr != ADJ_ERR_OK)
        return ierr;
    }
  }

  block_action_func(nb_variables, variables, dependencies, block.hermitian, block.coefficient, input, block.context, output );

  if (block.has_nonlinear_block)
    free(dependencies);

  return ADJ_ERR_OK;
}


int adj_evaluate_block_assembly(adj_adjointer* adjointer, adj_block block, adj_matrix *output, adj_vector* rhs)
{
  int i, ierr;
  void (*block_assembly_func)(int, adj_variable*, adj_vector*, int, adj_scalar, void*, adj_matrix*, adj_vector*) = NULL;
  adj_vector* dependencies = NULL;
  int nb_variables = 0;
  adj_variable* variables = NULL;

  ierr = adj_find_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, block. name, (void (**)(void)) &block_assembly_func);
  if (ierr != ADJ_ERR_OK)
    return ierr;

  /* We have the right callback, so let's call it already */ 
  if (block.has_nonlinear_block)
  {
    /* we need to set up the dependencies */
    nb_variables = block.nonlinear_block.ndepends;
    variables = block.nonlinear_block.depends;
    dependencies = (adj_vector*) malloc(nb_variables * sizeof(adj_vector));

    for (i = 0; i < nb_variables; i++)
    {
      ierr = adj_get_variable_value(adjointer, variables[i], &dependencies[i]);
      if (ierr != ADJ_ERR_OK)
        return ierr;
    }
  }

  block_assembly_func(nb_variables, variables, dependencies, block.hermitian, block.coefficient, block.context, output, rhs);

  if (block.has_nonlinear_block)
    free(dependencies);

  return ADJ_ERR_OK;
}


int adj_evaluate_nonlinear_derivative_action(adj_adjointer* adjointer, int nderivatives, adj_nonlinear_block_derivative* derivatives, adj_vector value, adj_vector* rhs)
{
  int ierr;
  int deriv;
  int rhs_allocated = 0;
  void (*nonlinear_derivative_action_func)(int nvar, adj_variable* variables, adj_vector* dependencies, adj_variable derivative, adj_vector contraction, int hermitian, adj_vector input, void* context, adj_vector* output);
  void (*nonlinear_action_func)(int nvar, adj_variable* variables, adj_vector* dependencies, adj_vector input, void* context, adj_vector* output) = NULL;

  /* As usual, check as much as we can at the start */
  strncpy(adj_error_msg, "Need a data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.vec_destroy == NULL)   return ADJ_ERR_NEED_CALLBACK;
  if (adjointer->callbacks.vec_axpy == NULL)      return ADJ_ERR_NEED_CALLBACK;
  if (adjointer->callbacks.vec_duplicate == NULL) return ADJ_ERR_NEED_CALLBACK;
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);
  assert(nderivatives > 0);

  /* Let's also check we have all of the variables available */
  for (deriv = 0; deriv < nderivatives; deriv++)
  {
    int i;
    for (i = 0; i < derivatives[deriv].nonlinear_block.ndepends; i++)
    {
      ierr = adj_has_variable_value(adjointer, derivatives[deriv].nonlinear_block.depends[i]);
      if (ierr != ADJ_ERR_OK) return ierr;
    }
  }

  for (deriv = 0; deriv < nderivatives; deriv++)
  {
    /* First, try to find a routine supplied by the user. */
    ierr = adj_find_operator_callback(adjointer, ADJ_NBLOCK_DERIVATIVE_ACTION_CB, derivatives[deriv].nonlinear_block.name, (void (**)(void)) &nonlinear_derivative_action_func);
    if (ierr == ADJ_ERR_OK)
    {
      ierr = adj_evaluate_nonlinear_derivative_action_supplied(adjointer, nonlinear_derivative_action_func, derivatives[deriv], value, rhs, &rhs_allocated);
      if (ierr != ADJ_ERR_OK) return ierr;
    }
    /* OK, if that didn't work, let's try ISP. */
    ierr = adj_find_operator_callback(adjointer, ADJ_NBLOCK_ACTION_CB, derivatives[deriv].nonlinear_block.name, (void (**)(void)) &nonlinear_action_func);
    if (ierr != ADJ_ERR_OK)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Could not find a nonlinear derivative action callback, nor a nonlinear action callback, for operator %s.", derivatives[deriv].nonlinear_block.name);
      return ierr;
    }
/*    ierr = adj_evaluate_nonlinear_derivative_action_isp(adjointer, nonlinear_action_func, derivatives[deriv], value, rhs, &rhs_allocated); */
    if (ierr != ADJ_ERR_OK) return ierr;
  }

  return ADJ_ERR_OK;
}

int adj_evaluate_nonlinear_derivative_action_supplied(adj_adjointer* adjointer, void (*nonlinear_derivative_action_func)(int nvar, adj_variable* variables, 
     adj_vector* dependencies, adj_variable derivative, adj_vector contraction, int hermitian, adj_vector input, void* context, adj_vector* output),
     adj_nonlinear_block_derivative derivative, adj_vector value, adj_vector* rhs, int* rhs_allocated)
{
  adj_vector* dependencies = NULL;
  adj_vector rhs_tmp;
  int i;
  int ierr;

  dependencies = (adj_vector*) malloc(derivative.nonlinear_block.ndepends * sizeof(adj_vector));
  for (i = 0; i < derivative.nonlinear_block.ndepends; i++)
  {
    ierr = adj_get_variable_value(adjointer, derivative.nonlinear_block.depends[i], &(dependencies[i]));
    assert(ierr == ADJ_ERR_OK); /* We checked for them earlier */
  }

  nonlinear_derivative_action_func(derivative.nonlinear_block.ndepends, derivative.nonlinear_block.depends, dependencies, derivative.variable,
                                   derivative.contraction, derivative.hermitian, value, derivative.nonlinear_block.context, &rhs_tmp);

  free(dependencies);

  if (!(*rhs_allocated))
  {
    adjointer->callbacks.vec_duplicate(rhs_tmp, rhs);
    *rhs_allocated = ADJ_TRUE;
  }

  adjointer->callbacks.vec_axpy(rhs, (adj_scalar) 1.0, rhs_tmp);
  adjointer->callbacks.vec_destroy(&rhs_tmp);

  return ADJ_ERR_OK;
}

/* int adj_evaluate_nonlinear_derivative_action_isp(adj_adjointer* adjointer, void (*nonlinear_action_func)(int nvar, adj_variable*
     variables, adj_vector* dependencies, adj_vector input, void* context, adj_vector* output), adj_nonlinear_block_derivative derivative, 
     adj_vector value, adj_vector* rhs, int* rhs_allocated)
{
} */

int adj_evaluate_functional(adj_adjointer* adjointer, adj_variable variable, char* functional, adj_vector* output)
{
  int i, ierr;
  void (*functional_derivative_func)(adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar start_time, adj_scalar end_time, adj_vector* output) = NULL;
  adj_vector* dependencies = NULL;
  int nb_variables = 0;
  adj_variable* variables = NULL;
  adj_functional_data* functional_data_ptr = NULL;
  adj_scalar start_time, end_time;

  ierr = adj_find_functional_derivative_callback(adjointer, functional, &functional_derivative_func);
  if (ierr != ADJ_ERR_OK)
    return ierr;

  if (adjointer->ntimesteps <= variable.timestep) 
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "No data is associated with this timestep %d.", variable.timestep);
    return ADJ_ERR_INVALID_INPUTS;
  }

  functional_data_ptr = adjointer->timestep_data[variable.timestep].functional_data_start;
  while (functional_data_ptr != NULL)
  {
    if (strncmp(functional_data_ptr->name, functional, ADJ_NAME_LEN) == 0)
    {
      nb_variables = functional_data_ptr->ndepends;
      variables = functional_data_ptr->dependencies;
      dependencies = (adj_vector*) malloc(nb_variables * sizeof(adj_vector));
      /* we need to set up the dependencies */
      for (i = 0; i < nb_variables; i++)
      {
        ierr = adj_get_variable_value(adjointer, variables[i], &dependencies[i]);
        if (ierr != ADJ_ERR_OK)
          return ierr;
      }
      break;
    }
    functional_data_ptr = functional_data_ptr->next;
  }

  start_time = adjointer->timestep_data[variable.timestep].start_time;
  end_time = adjointer->timestep_data[variable.timestep].end_time;

  /* We have the right callback, so let's call it already */ 
  functional_derivative_func(variable, nb_variables, variables, dependencies, functional, start_time, end_time, output);

  if (dependencies != NULL) free(dependencies);

  return ADJ_ERR_OK;
}
