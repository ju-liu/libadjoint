#include "libadjoint/adj_evaluation.h"

int adj_evaluate_block_action(adj_adjointer* adjointer, adj_block block, adj_vector input, adj_vector* output)
{
  int i, ierr;
  void (*block_action_func)(int, adj_variable*, adj_vector*, int, adj_vector, void*, adj_vector*)=NULL;
  adj_vector* dependencies = NULL;
  int nb_variables=0;
  adj_variable* variables=NULL;

  ierr = adj_find_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, block.name, (void (**)(void)) &block_action_func);
  if (ierr!=0)
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
  block_action_func(nb_variables, variables, dependencies, block.hermitian, input, block.context, output );
  return 0;
}


int adj_evaluate_block_assembly(adj_adjointer* adjointer, adj_block block, adj_matrix *output)
{
  int i, ierr;
  void (*block_assembly_func)(int, adj_variable*, adj_vector*, int, void*, adj_matrix*)=NULL;
  adj_vector* dependencies = NULL;
  int nb_variables=0;
  adj_variable* variables=NULL;

  ierr = adj_find_operator_callback(adjointer, ADJ_BLOCK_ASSEMBLY_CB, block. name, (void (**)(void)) &block_assembly_func);
  if (ierr!=0)
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
  block_assembly_func(nb_variables, variables, dependencies, block.hermitian, block.context, output );
  return 0;
}


/* int adj_evaluate_nonlinear_derivative_action(adj_adjointer* adjointer, adj_nonlinear_block_derivative* derivatives, adj_vector value, adj_vector rhs)
{
 
} */
