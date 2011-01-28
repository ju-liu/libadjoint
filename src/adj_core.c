#include "libadjoint/adj_core.h"

int adj_get_adjoint_equation(adj_adjointer* adjointer, int equation, int functional, adj_matrix* lhs, adj_vector* rhs, adj_variable* adj_var)
{
  int ierr;
  adj_equation fwd_eqn;
  adj_variable fwd_var;
  adj_variable_data* adj_data;
  adj_variable_data* fwd_data;
  int i;
  int j;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING)
  {
    strncpy(adj_error_msg, "You have asked for an adjoint equation, but the adjointer has been deactivated.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  if (equation < 0 || equation > adjointer->nequations)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Invalid equation number %d.", equation);
    return ADJ_ERR_INVALID_INPUTS;
  }

  fwd_eqn = adjointer->equations[equation];
  fwd_var = fwd_eqn.variable;

  /* Create the associated adjoint variable */
  ierr = adj_create_variable(fwd_var.name, fwd_var.timestep, fwd_var.iteration, fwd_var.auxiliary, adj_var);
  if (ierr != ADJ_ERR_OK) return ierr;
  adj_var->type = ADJ_ADJOINT;
  adj_var->functional = functional;

  /* Add an entry in the hash table for this variable */
  ierr = adj_find_variable_data(&(adjointer->varhash), adj_var, &adj_data);
  if (ierr == ADJ_ERR_HASH_FAILED)
  {
    /* It might not fail, if we have tried to fetch this equation already */
    ierr = adj_add_new_hash_entry(adjointer, adj_var, &adj_data);
    if (ierr != ADJ_ERR_OK) return ierr;
  }

  /* Now let's fill in its data */
  adj_data->equation = -1; /* it never has a forward equation */
  /* And fill in its .adjoint_equations */
  /* The adjoint equations this variable is necessary for are:
     * The (adjoint equation) of (the target) of (each block) in (the forward equation) associated with (the adjoint equation we're fetching)
     * The (adjoint equation) of (the dependencies) of (each block) in (the forward equation) associated with (the adjoint equation we're fetching)
   Do you see why working that out gave me an almighty headache? */
  for (i = 0; i < fwd_eqn.nblocks; i++)
  {
    adj_variable_data* block_target_data;
    ierr = adj_find_variable_data(&(adjointer->varhash), &(fwd_eqn.targets[i]), &block_target_data);
    if (ierr != ADJ_ERR_OK) return ierr;
    adj_append_unique(&(adj_data->adjoint_equations), &(adj_data->nadjoint_equations), block_target_data->equation);

    for (j = 0; j < fwd_eqn.blocks[i].nonlinear_block.ndepends; j++)
    {
      adj_variable_data* j_data;
      ierr = adj_find_variable_data(&(adjointer->varhash), &(fwd_eqn.blocks[i].nonlinear_block.depends[j]), &j_data);
      if (ierr != ADJ_ERR_OK) return ierr;
      adj_append_unique(&(adj_data->adjoint_equations), &(adj_data->nadjoint_equations), j_data->equation);
    }
  }

  /* --------------------------------------------------------------------------
   * Computation of A* terms                                                  |
   * -------------------------------------------------------------------------- */

  ierr = adj_find_variable_data(&(adjointer->varhash), &fwd_var, &fwd_data);
  if (ierr != ADJ_ERR_OK) return ierr;

  /* fwd_data->targeting_equations what forward equations have nonzero blocks in the column of A associated with fwd_var. */
  /* That column of A becomes the current row of A* we want to now compute. */

  /* First we find the diagonal entry and assemble that, to compute the A* of lhs. */
  /* Find the block in fwd_eqn that targets fwd_var, and assemble that (hermitianed) */
  {
    adj_block block;
    for (i = 0; i < fwd_eqn.nblocks; i++)
    {
      if (adj_variable_equal(&(fwd_eqn.targets[i]), &fwd_var, 1))
      {
        /* this is the right block */
        block = fwd_eqn.blocks[i];
        break;
      }
    }
    block.hermitian = 1;
    ierr = adj_evaluate_block_assembly(adjointer, block, lhs);
    if (ierr != ADJ_ERR_OK) return ierr;
  }

  /* Great! Now let's assemble the RHS contributions of A*. */

  return ADJ_ERR_OK;
}
