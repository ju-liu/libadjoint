#include "libadjoint/adj_core.h"

int adj_get_adjoint_equation(adj_adjointer* adjointer, int equation, char* functional, adj_matrix* lhs, adj_vector* rhs, adj_variable* adj_var)
{
  int ierr;
  adj_equation fwd_eqn;
  adj_variable fwd_var;
  adj_variable_data* adj_data;
  adj_variable_data* fwd_data;
  int i;
  int j;
  void (*functional_derivative_func)(adj_adjointer* adjointer, adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output) = NULL;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING)
  {
    strncpy(adj_error_msg, "You have asked for an adjoint equation, but the adjointer has been deactivated.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  if (equation < 0 || equation >= adjointer->nequations)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Invalid equation number %d.", equation);
    return ADJ_ERR_INVALID_INPUTS;
  }

  strncpy(adj_error_msg, "Need a data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.vec_destroy == NULL) return ADJ_ERR_NEED_CALLBACK;
  if (adjointer->callbacks.vec_axpy == NULL)    return ADJ_ERR_NEED_CALLBACK;
  if (adjointer->callbacks.mat_axpy == NULL)    return ADJ_ERR_NEED_CALLBACK;
  if (adjointer->callbacks.mat_destroy == NULL) return ADJ_ERR_NEED_CALLBACK;
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  ierr = adj_find_functional_derivative_callback(adjointer, functional, &functional_derivative_func);
  if (ierr != ADJ_ERR_OK) return ierr;

  fwd_eqn = adjointer->equations[equation];
  fwd_var = fwd_eqn.variable;

  ierr = adj_find_variable_data(&(adjointer->varhash), &fwd_var, &fwd_data);
  assert(ierr == ADJ_ERR_OK);

  /* Check that we have all the adjoint values we need, before we start allocating stuff */
  for (i = 0; i < fwd_data->ntargeting_equations; i++)
  {
    adj_equation other_fwd_eqn;
    adj_variable other_adj_var;

    if (fwd_data->targeting_equations[i] == equation) continue; /* that term goes in the lhs, and we've already taken care of it */
    other_fwd_eqn = adjointer->equations[fwd_data->targeting_equations[i]];

    /* Find the adjoint variable we want this to multiply */
    other_adj_var = other_fwd_eqn.variable; other_adj_var.type = ADJ_ADJOINT; strncpy(other_adj_var.functional, functional, ADJ_NAME_LEN);
    /* and now get its value */
    ierr = adj_has_variable_value(adjointer, other_adj_var);
    if (ierr != ADJ_ERR_OK)
    {
      char buf[255];
      adj_variable_str(other_adj_var, buf, 255);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for variable %s, but don't have one.", buf);
      return ADJ_ERR_NEED_VALUE;
    }
  }

  /* Create the associated adjoint variable */
  ierr = adj_create_variable(fwd_var.name, fwd_var.timestep, fwd_var.iteration, fwd_var.auxiliary, adj_var);
  if (ierr != ADJ_ERR_OK) return ierr;
  adj_var->type = ADJ_ADJOINT;
  strncpy(adj_var->functional, functional, ADJ_NAME_LEN);

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
    block.hermitian = !block.hermitian;
    ierr = adj_evaluate_block_assembly(adjointer, block, lhs, rhs);
    if (ierr != ADJ_ERR_OK) return ierr;
  }

  /* Great! Now let's assemble the RHS contributions of A*. */

  /* Now loop through the off-diagonal blocks of A*. */
  for (i = 0; i < fwd_data->ntargeting_equations; i++)
  {
    adj_equation other_fwd_eqn;
    adj_block block;
    adj_variable other_adj_var;
    adj_vector adj_value;
    adj_vector rhs_tmp;

    if (fwd_data->targeting_equations[i] == equation) continue; /* that term goes in the lhs, and we've already taken care of it */
    other_fwd_eqn = adjointer->equations[fwd_data->targeting_equations[i]];

    for (j = 0; j < other_fwd_eqn.nblocks; j++)
    {
      if (adj_variable_equal(&fwd_var, &(other_fwd_eqn.targets[j]), 1))
      {
        block = other_fwd_eqn.blocks[j];
        break;
      }
    }

    /* OK. Now we've found the right block ... */
    block.hermitian = !block.hermitian;

    /* Find the adjoint variable we want this to multiply */
    other_adj_var = other_fwd_eqn.variable; other_adj_var.type = ADJ_ADJOINT; strncpy(other_adj_var.functional, functional, ADJ_NAME_LEN);
    /* and now get its value */
    ierr = adj_get_variable_value(adjointer, other_adj_var, &adj_value);
    assert(ierr == ADJ_ERR_OK); /* we should have them all, we checked for them earlier */

    ierr = adj_evaluate_block_action(adjointer, block, adj_value, &rhs_tmp);
    if (ierr != ADJ_ERR_OK) return ierr;
    adjointer->callbacks.vec_axpy(rhs, (adj_scalar)-1.0, rhs_tmp);
    adjointer->callbacks.vec_destroy(&rhs_tmp);
  }

  /* Now add dJ/du to the rhs */
  {
    adj_vector rhs_tmp;
    int has_djdu;
    ierr = adj_evaluate_functional_derivative(adjointer, fwd_var, functional, &rhs_tmp, &has_djdu);
    if (ierr != ADJ_ERR_OK) return ierr;
    if (has_djdu)
    {
      adjointer->callbacks.vec_axpy(rhs, (adj_scalar)1.0, rhs_tmp);
      adjointer->callbacks.vec_destroy(&rhs_tmp);
    }
  }

  return ADJ_ERR_OK;
}

int adj_get_forward_equation(adj_adjointer* adjointer, int equation, adj_matrix* lhs, adj_vector* rhs, adj_variable* fwd_var)
{
  int ierr;
  adj_equation fwd_eqn;
  adj_variable_data* fwd_data;
  adj_vector rhs_tmp;
  int i;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING)
  {
    strncpy(adj_error_msg, "You have asked for an forward equation, but the adjointer has been deactivated.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  if (equation < 0 || equation >= adjointer->nequations)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Invalid equation number %d.", equation);
    return ADJ_ERR_INVALID_INPUTS;
  }

  strncpy(adj_error_msg, "Need a data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.vec_destroy == NULL) return ADJ_ERR_NEED_CALLBACK;
  if (adjointer->callbacks.vec_axpy == NULL)    return ADJ_ERR_NEED_CALLBACK;
  if (adjointer->callbacks.mat_axpy == NULL)    return ADJ_ERR_NEED_CALLBACK;
  if (adjointer->callbacks.mat_destroy == NULL) return ADJ_ERR_NEED_CALLBACK;
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  fwd_eqn = adjointer->equations[equation];
  *fwd_var = fwd_eqn.variable;

  ierr = adj_find_variable_data(&(adjointer->varhash), fwd_var, &fwd_data);
  assert(ierr == ADJ_ERR_OK);

  /* Check that we have all the forward values we need, before we start allocating stuff */
  for (i = 0; i < fwd_eqn.nblocks; i++)
  {
    adj_variable other_adj_var;

    /* Get the forward variable we want this to multiply */
    other_adj_var = fwd_eqn.targets[i];
    if (adj_variable_equal(fwd_var, &other_adj_var, 1)) continue; /* that term goes in the lhs */
    /* and now check it has a value */
    ierr = adj_has_variable_value(adjointer, other_adj_var);
    if (ierr != ADJ_ERR_OK)
    {
      char buf[255];
      adj_variable_str(other_adj_var, buf, 255);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for variable %s, but don't have one.", buf);
      return ADJ_ERR_NEED_VALUE;
    }
  }

  /* --------------------------------------------------------------------------
   * Computation of A terms                                                  |
   * -------------------------------------------------------------------------- */

  /* fwd_data->targeting_equations what forward equations have nonzero blocks in the column of A associated with fwd_var. */

  {
    adj_block block;
    for (i = 0; i < fwd_eqn.nblocks; i++)
    {
      if (adj_variable_equal(&(fwd_eqn.targets[i]), fwd_var, 1))
      {
        /* this is the right block */
        block = fwd_eqn.blocks[i];
        break;
      }
    }
    ierr = adj_evaluate_block_assembly(adjointer, block, lhs, rhs);
    if (ierr != ADJ_ERR_OK) return ierr;
  }

  /* Great! Now let's assemble the RHS contributions of A. */

  /* Now loop through the off-diagonal blocks of A. */
  for (i = 0; i < fwd_eqn.nblocks; i++)
  {
    adj_block block;
    adj_variable other_var;
    adj_vector value;

    /* Get the forward variable we want this block to multiply with */
    other_var = fwd_eqn.targets[i];

    /* Ignore the diagonal block */
    if (adj_variable_equal(&other_var, fwd_var, 1))
      continue;

    block = fwd_eqn.blocks[i];

    /* and now get its value */
    ierr = adj_get_variable_value(adjointer, other_var, &value);
    assert(ierr == ADJ_ERR_OK); /* we should have them all, we checked for them earlier */

    ierr = adj_evaluate_block_action(adjointer, block, value, &rhs_tmp);
    if (ierr != ADJ_ERR_OK) return ierr;
    adjointer->callbacks.vec_axpy(rhs, (adj_scalar)-1.0, rhs_tmp);
    adjointer->callbacks.vec_destroy(&rhs_tmp);
  }
  return ADJ_ERR_OK;

}

