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


/* int adj_evaluate_nonlinear_derivative_action(adj_adjointer* adjointer, adj_nonlinear_block_derivative* derivatives, adj_vector value, adj_vector rhs)
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

  if (nb_variables > 0)
    free(dependencies);

  return ADJ_ERR_OK;
}
