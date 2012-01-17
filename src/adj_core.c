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
  void (*functional_derivative_func)(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output) = NULL;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING)
  {
    strncpy(adj_error_msg, "You have asked for an adjoint equation, but the adjointer has been deactivated.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (equation < 0 || equation >= adjointer->nequations)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Invalid equation number %d.", equation);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (adjointer->callbacks.vec_destroy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_VEC_DESTROY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.vec_axpy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_VEC_AXPY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.mat_axpy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_MAT_AXPY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.mat_destroy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_MAT_DESTROY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }

  ierr = adj_find_functional_derivative_callback(adjointer, functional, &functional_derivative_func);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

  fwd_eqn = adjointer->equations[equation];
  fwd_var = fwd_eqn.variable;

  ierr = adj_find_variable_data(&(adjointer->varhash), &fwd_var, &fwd_data);
  assert(ierr == ADJ_OK);

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
    if (ierr != ADJ_OK)
    {
      char buf[255];
      adj_variable_str(other_adj_var, buf, 255);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for variable %s, but don't have one.", buf);
      return adj_chkierr_auto(ADJ_ERR_NEED_VALUE);
    }
  }

  /* Create the associated adjoint variable */
  ierr = adj_create_variable(fwd_var.name, fwd_var.timestep, fwd_var.iteration, fwd_var.auxiliary, adj_var);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
  adj_var->type = ADJ_ADJOINT;
  strncpy(adj_var->functional, functional, ADJ_NAME_LEN);

  /* Add an entry in the hash table for this variable */
  ierr = adj_find_variable_data(&(adjointer->varhash), adj_var, &adj_data);
  if (ierr == ADJ_ERR_HASH_FAILED)
  {
    /* It might not fail, if we have tried to fetch this equation already */
    ierr = adj_add_new_hash_entry(adjointer, adj_var, &adj_data);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
  }

  /* Now let's fill in its data */
  adj_data->equation = -1; /* it never has a forward equation */
  /* And fill in its .adjoint_equations */
  /* The adjoint equations this variable is necessary for are:
     * The (adjoint equation) of (the target) of (each block) in (the forward equation) associated with (the adjoint equation we're fetching)
     * The (adjoint equation) of (the dependencies) of (each block) in (the forward equation) associated with (the adjoint equation we're fetching)
     * The (adjoint equation) of (the dependencies) of (the right-hand-side) of (the forward equation) associated with (the adjoint equation we're fetching)
   Do you see why working that out gave me an almighty headache? */
  for (i = 0; i < fwd_eqn.nblocks; i++)
  {
    /* A* terms */
    adj_variable_data* block_target_data;
    ierr = adj_find_variable_data(&(adjointer->varhash), &(fwd_eqn.targets[i]), &block_target_data);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adj_append_unique(&(adj_data->adjoint_equations), &(adj_data->nadjoint_equations), block_target_data->equation);

    /* G* terms */
    for (j = 0; j < fwd_eqn.blocks[i].nonlinear_block.ndepends; j++)
    {
      adj_variable_data* j_data;
      ierr = adj_find_variable_data(&(adjointer->varhash), &(fwd_eqn.blocks[i].nonlinear_block.depends[j]), &j_data);
      if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      adj_append_unique(&(adj_data->adjoint_equations), &(adj_data->nadjoint_equations), j_data->equation);
    }
  }
  /* R* terms */
  for (i = 0; i < fwd_eqn.nrhsdeps; i++)
  {
    adj_variable_data* rhs_dep_data;
    ierr = adj_find_variable_data(&(adjointer->varhash), &(fwd_eqn.rhsdeps[i]), &rhs_dep_data);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adj_append_unique(&(adj_data->adjoint_equations), &(adj_data->nadjoint_equations), rhs_dep_data->equation);
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
    int blockcount = 0;
    for (i = 0; i < fwd_eqn.nblocks; i++)
    {

      if (adj_variable_equal(&(fwd_eqn.targets[i]), &fwd_var, 1))
      {
        /* this is the right block */
        block = fwd_eqn.blocks[i];
        block.hermitian = !block.hermitian;
        blockcount++;
        if (blockcount == 1) /* the first one we've found */
        {
          ierr = adj_evaluate_block_assembly(adjointer, block, lhs, rhs);
          if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        }
        else
        {
          adj_matrix lhs_tmp;
          adj_vector rhs_tmp;
          ierr = adj_evaluate_block_assembly(adjointer, block, &lhs_tmp, &rhs_tmp);
          if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
          adjointer->callbacks.vec_destroy(&rhs_tmp); /* we already have rhs from the first block assembly */
          adjointer->callbacks.mat_axpy(lhs, (adj_scalar) 1.0, lhs_tmp); /* add lhs_tmp to lhs */
          adjointer->callbacks.mat_destroy(&lhs_tmp);
        }
      }
    }
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

        /* OK. Now we've found the right block ... */
        block.hermitian = !block.hermitian;

        /* Find the adjoint variable we want this to multiply */
        other_adj_var = other_fwd_eqn.variable; other_adj_var.type = ADJ_ADJOINT; strncpy(other_adj_var.functional, functional, ADJ_NAME_LEN);
        /* and now get its value */
        ierr = adj_get_variable_value(adjointer, other_adj_var, &adj_value);
        assert(ierr == ADJ_OK); /* we should have them all, we checked for them earlier */

        ierr = adj_evaluate_block_action(adjointer, block, adj_value, &rhs_tmp);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        adjointer->callbacks.vec_axpy(rhs, (adj_scalar)-1.0, rhs_tmp);
        adjointer->callbacks.vec_destroy(&rhs_tmp);
      }
    }

  }

  /* --------------------------------------------------------------------------
   * Computation of G* terms                                                  |
   * -------------------------------------------------------------------------- */

  /* We need to loop through the equations that depend on fwd_var; each one of those will produce
     a term in this row of G*. */
  {
    for (i = 0; i < fwd_data->ndepending_equations; i++)
    {
      /* We compute the block of G*, either assembling it, or computing its action, as appropriate. */
      int ndepending_eqn;
      adj_equation depending_eqn;
      int nderivs; /* these two are the raw derivatives to compute */
      adj_nonlinear_block_derivative* derivs;

      int nnew_derivs; /* and these two are after derivative simplification */
      adj_nonlinear_block_derivative* new_derivs;
      int l, k;

      ndepending_eqn = fwd_data->depending_equations[i];
      depending_eqn = adjointer->equations[ndepending_eqn];

      /* Now, we loop through the blocks of this forward equation, finding the ones that depend on the fwd variable.
         First, let's find out how many blocks depend on this variable; then we'll malloc that many 
         adj_nonlinear_block_derivatives, and then we'll go about filling them in. */
      nderivs = 0;
      for (j = 0; j < depending_eqn.nblocks; j++)
      {
        if (depending_eqn.blocks[j].has_nonlinear_block)
        {
          for (k = 0; k < depending_eqn.blocks[j].nonlinear_block.ndepends; k++)
          {
            if (adj_variable_equal(&fwd_var, &depending_eqn.blocks[j].nonlinear_block.depends[k], 1))
            {
              nderivs++;
            }
          }
        }
      }

      /* Now that we have ndepending_blocks, let's use it */
      derivs = (adj_nonlinear_block_derivative*) malloc(nderivs * sizeof(adj_nonlinear_block_derivative));
      ADJ_CHKMALLOC(derivs);
      l = 0;
      for (j = 0; j < depending_eqn.nblocks; j++)
      {
        if (depending_eqn.blocks[j].has_nonlinear_block)
        {
          for (k = 0; k < depending_eqn.blocks[j].nonlinear_block.ndepends; k++)
          {
            if (adj_variable_equal(&fwd_var, &depending_eqn.blocks[j].nonlinear_block.depends[k], 1))
            {
              adj_vector target;
              ierr = adj_get_variable_value(adjointer, depending_eqn.targets[j], &target);
              if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

              ierr = adj_create_nonlinear_block_derivative(adjointer, depending_eqn.blocks[j].nonlinear_block, depending_eqn.blocks[j].coefficient, fwd_var, target, !depending_eqn.blocks[j].hermitian, &derivs[l]);
              if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
              l++;
            }
          }
        }
      }

      /* OK, Here's where we do our simplifications; this can be a significant optimisation */
      /* .......................................................................................... */
      ierr = adj_simplify_derivatives(adjointer, nderivs, derivs, &nnew_derivs, &new_derivs);
      for (l = 0; l < nderivs; l++)
      {
        ierr = adj_destroy_nonlinear_block_derivative(adjointer, &derivs[l]);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      }
      free(derivs);

      /* Now, we go and evaluate each one of the derivatives, assembling or acting as appropriate */
      if (ndepending_eqn == equation)
      {
        /* This G-block is on the diagonal, and so we must assemble it .. later */
        snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Sorry, we can't handle G-blocks on the diagonal (yet).");
        return adj_chkierr_auto(ADJ_ERR_NOT_IMPLEMENTED);
      }
      else
      {
        /* This G-block is NOT on the diagonal, so we only need its action */
        adj_variable adj_associated;
        adj_vector adj_value;

        adj_associated = depending_eqn.variable;
        adj_associated.type = ADJ_ADJOINT;
        strncpy(adj_associated.functional, functional, ADJ_NAME_LEN);
        ierr = adj_get_variable_value(adjointer, adj_associated, &adj_value);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        /* And now we are ready */
        ierr = adj_evaluate_nonlinear_derivative_action(adjointer, nnew_derivs, new_derivs, adj_value, rhs);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      }

      for (l = 0; l < nnew_derivs; l++)
      {
        ierr = adj_destroy_nonlinear_block_derivative(adjointer, &new_derivs[l]);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      }
      free(new_derivs);
    }
  }

  /* --------------------------------------------------------------------------
   * Computation of R* terms                                                  |
   * -------------------------------------------------------------------------- */
  {
    for (i = 0; i < fwd_data->nrhs_equations; i++)
    {
      int rhs_equation = fwd_data->rhs_equations[i];
      /* Does this R* contribute to the adjoint matrix ... */
      if (adj_variable_equal(&(adjointer->equations[rhs_equation].variable), &fwd_var, 1))
      {
        adj_matrix rstar;
        ierr = adj_evaluate_rhs_deriv_assembly(adjointer, adjointer->equations[rhs_equation], ADJ_TRUE, &rstar);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        adjointer->callbacks.mat_axpy(lhs, (adj_scalar) -1.0, rstar); /* Subtract the R* contribution from the adjoint lhs */
        adjointer->callbacks.mat_destroy(&rstar);
      }
      /* ... or to the right-hand side of the adjoint system? */
      else
      {
        /* Get the adj_equation associated with this dependency, so we can pull out the relevant rhs_deriv_action callback */
        adj_vector deriv_action;
        adj_variable contraction_var;
        adj_vector contraction;
        int has_output;

        has_output = -666;

        contraction_var = adjointer->equations[rhs_equation].variable; contraction_var.type = ADJ_ADJOINT; strncpy(contraction_var.functional, functional, ADJ_NAME_LEN);
        ierr = adj_get_variable_value(adjointer, contraction_var, &contraction);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        ierr = adj_evaluate_rhs_deriv_action(adjointer, adjointer->equations[rhs_equation], fwd_var, contraction, ADJ_TRUE, &deriv_action, &has_output);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        if (has_output)
        {
          /* Now that we have the contribution, we need to add it to the adjoint right hand side */
          adjointer->callbacks.vec_axpy(rhs, (adj_scalar)1.0, deriv_action);
          adjointer->callbacks.vec_destroy(&deriv_action);
        }
      }
    }
  }

  /* Now add dJ/du to the rhs */
  {
    adj_vector rhs_tmp;
    int has_djdu;
    ierr = adj_evaluate_functional_derivative(adjointer, fwd_var, functional, &rhs_tmp, &has_djdu);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    if (has_djdu)
    {
      adjointer->callbacks.vec_axpy(rhs, (adj_scalar)1.0, rhs_tmp);
      adjointer->callbacks.vec_destroy(&rhs_tmp);
    }
  }

  return ADJ_OK;
}

int adj_get_adjoint_solution(adj_adjointer* adjointer, int equation, char* functional, adj_vector* soln, adj_variable* adj_var)
{
  int ierr;
  adj_matrix lhs;
  adj_vector rhs;
  int cs;

  /* Check for the required callbacks */ 
  if (adjointer->callbacks.solve == NULL)
  {   
    strncpy(adj_error_msg, "Need the solve data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }

  /* If using revolve, we might have to restore from a checkpoint before the adjoint equation can be solved */
  ierr = adj_get_checkpoint_strategy(adjointer, &cs);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

  /* Online revolve needs to be told the total number of timesteps
   * before solving the first adjoint equation */
  if (cs == ADJ_CHECKPOINT_REVOLVE_ONLINE && equation == adjointer->nequations-1)
    {
      adjointer->revolve_data.steps = adjointer->nequations;
      revolve_turn(adjointer->revolve_data.revolve, adjointer->revolve_data.steps);
    }

  if ((cs == ADJ_CHECKPOINT_REVOLVE_OFFLINE) || (cs == ADJ_CHECKPOINT_REVOLVE_MULTISTAGE) || (cs == ADJ_CHECKPOINT_REVOLVE_ONLINE))
  {
    /* Recompute any variables that is required for solving this adjoint equation */
    ierr = adj_revolve_to_adjoint_equation(adjointer, equation);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
  }

  if (adjointer->revolve_data.verbose)
     printf("Revolve: Solving adjoint equation %i.\n", equation);

  /* At this point, all the dependencies are available to assemble the adjoint equation */
  ierr = adj_get_adjoint_equation(adjointer, equation, functional, &lhs, &rhs, adj_var);
  if (ierr != ADJ_OK)
    return adj_chkierr_auto(ierr);

  /* Solve the linear system */
  adjointer->callbacks.solve(*adj_var, lhs, rhs, soln); 
  adjointer->callbacks.vec_destroy(&rhs);
  adjointer->callbacks.mat_destroy(&lhs);

  /* We can now safely un-checkoint this equation and its associated forward variable */
  if ((cs == ADJ_CHECKPOINT_REVOLVE_OFFLINE) || (cs == ADJ_CHECKPOINT_REVOLVE_MULTISTAGE) || (cs == ADJ_CHECKPOINT_REVOLVE_ONLINE))
  {
    adj_variable_data* data_ptr;

    ierr = adj_find_variable_data(&(adjointer->varhash), &adjointer->equations[equation].variable, &data_ptr);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

    if (adjointer->revolve_data.verbose)
      if ((adjointer->equations[equation].disk_checkpoint == ADJ_TRUE) ||
          (adjointer->equations[equation].memory_checkpoint == ADJ_TRUE))
        printf("Revolve: Delete checkpoint equation %i.\n", equation);

    adjointer->equations[equation].disk_checkpoint=ADJ_FALSE;
    adjointer->equations[equation].memory_checkpoint=ADJ_FALSE;
    data_ptr->storage.storage_memory_is_checkpoint=ADJ_FALSE;
    data_ptr->storage.storage_disk_is_checkpoint=ADJ_FALSE;
  }

  return ADJ_OK;
}

int adj_revolve_to_adjoint_equation(adj_adjointer* adjointer, int equation)
{
  int ierr, cs;
  int capo, oldcapo;
  int start_eqn, end_eqn;
  int loop = ADJ_TRUE;

  ierr = adj_get_checkpoint_strategy(adjointer, &cs);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

  while(loop)
  {
    switch (adjointer->revolve_data.current_action)
    {
      case CACTION_ADVANCE:
        oldcapo = revolve_getoldcapo(adjointer->revolve_data.revolve);
        capo = revolve_getcapo(adjointer->revolve_data.revolve);

        /* If revolve advances to a timestep larger than ntimeteps-1,
         * then the user claimed in the revolve settings to solve for more timesteps then we actually did.
         */
        if (capo > adjointer->ntimesteps-1)
        {
          snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Revolve expected %i timesteps but only %i were annotated.", capo+1, adjointer->ntimesteps);
          return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
        }

        ierr = adj_timestep_start_equation(adjointer, oldcapo, &start_eqn);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        ierr = adj_timestep_end_equation(adjointer, capo-1, &end_eqn);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        if (adjointer->revolve_data.verbose)
          printf("====== Revolve: Replay from equation %i (first equation of timestep %i) to equation %i (last equation of timestep %i) =======\n", start_eqn, oldcapo, end_eqn, capo-1);

        ierr = adj_replay_forward_equations(adjointer, start_eqn, end_eqn, ADJ_FALSE);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        adjointer->revolve_data.current_timestep = capo;
        adjointer->revolve_data.current_action = revolve(adjointer->revolve_data.revolve);
        assert((adjointer->revolve_data.current_action == CACTION_TAKESHOT) ||
               (adjointer->revolve_data.current_action == CACTION_YOUTURN) ||
               (adjointer->revolve_data.current_action == CACTION_FIRSTRUN));

        break;

      case CACTION_TAKESHOT:
        /* Record all required forward variables for the first equation of the current timestep */
        ierr = adj_timestep_start_equation(adjointer, adjointer->revolve_data.current_timestep, &start_eqn);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        if (adjointer->revolve_data.verbose)
          printf("Revolve: Create checkpoint of equation %i (first equation of timestep %i).\n", start_eqn, adjointer->revolve_data.current_timestep);

        /* in a multistage setting, we have to ask revolve where to store the checkpoint */
        if ((cs == ADJ_CHECKPOINT_REVOLVE_MULTISTAGE) && (revolve_getwhere(adjointer->revolve_data.revolve) == 1))
          ierr = adj_checkpoint_equation(adjointer, start_eqn, ADJ_CHECKPOINT_STORAGE_MEMORY);
         else
          ierr = adj_checkpoint_equation(adjointer, start_eqn, ADJ_CHECKPOINT_STORAGE_DISK);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        adjointer->revolve_data.current_action = revolve(adjointer->revolve_data.revolve);
        break;

      case CACTION_FIRSTRUN:
        /* Check that the forward simulation was run to the last timestep */
        if (adjointer->revolve_data.current_timestep != adjointer->revolve_data.steps-1)
        {
          snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You asked for an adjoint solution after solving %i forward timestep, but you told revolve that you are going to solve %i forward timesteps.", adjointer->revolve_data.current_timestep, adjointer->revolve_data.steps);
          return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
        }

        /* Replay the last equation (if not already recorded)
         * This is done by leaving out the break here,
         * which executes the code for the CACTION_YOUTURN
         */

      case CACTION_YOUTURN:
        ierr = adj_timestep_start_equation(adjointer, adjointer->revolve_data.current_timestep, &start_eqn);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        ierr = adj_timestep_end_equation(adjointer, adjointer->revolve_data.current_timestep, &end_eqn);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        /* When replaying the timestep of the current adjoint equation, the replay records all variables of that timestep. */
        /* For that reason, we need to execute the replay only if we are about to solve the last equation of a timestep */
        if (equation == end_eqn)
        {
          if (adjointer->revolve_data.verbose)
            printf("====== Revolve: Replay from equation %i (first equation of timestep %i) to equation %i (last equation of timestep %i). ======\n", start_eqn, adjointer->revolve_data.current_timestep, end_eqn, adjointer->revolve_data.current_timestep);

          /* While replaying, we want to store the solved variables as checkpoints to ensure that we have all variables available for the upcoming adjoint solve */
          ierr = adj_replay_forward_equations(adjointer, start_eqn, end_eqn, ADJ_TRUE);
          if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        }

        loop=ADJ_FALSE;
        break;

      case CACTION_RESTORE:
        adjointer->revolve_data.current_timestep = revolve_getcapo(adjointer->revolve_data.revolve);
        adjointer->revolve_data.current_action = revolve(adjointer->revolve_data.revolve);
        break;

      case CACTION_ERROR:
        snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "An internal error occured: Irregular termination of revolve.");
        return adj_chkierr_auto(ADJ_ERR_REVOLVE_ERROR);

      default:
        snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "An internal error occured: The adjointer and revolve are out of sync.");
        return adj_chkierr_auto(ADJ_ERR_REVOLVE_ERROR);        
    }

  }

  /* Check that the timesteps are in sync */
  if (adjointer->revolve_data.current_timestep != adjointer->equations[equation].variable.timestep)
  {
    strncpy(adj_error_msg, "You must loop over the adjoint equation chronologically backwards in time, but revolve thinks you do not.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  /* If this function was called just before solving the last adjoint equation of the current timestep, then we ask revolve what to do next */
  if (equation == 0 || adjointer->revolve_data.current_timestep != adjointer->equations[equation-1].variable.timestep)
  {
    adjointer->revolve_data.current_action = revolve(adjointer->revolve_data.revolve);
  }
  return ADJ_OK;
}


int adj_get_forward_equation(adj_adjointer* adjointer, int equation, adj_matrix* lhs, adj_vector* rhs, adj_variable* fwd_var)
{
  int ierr;
  adj_equation fwd_eqn;
  adj_variable_data* fwd_data;
  adj_vector rhs_tmp;
  int i, j;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING)
  {
    strncpy(adj_error_msg, "You have asked for an forward equation, but the adjointer has been deactivated.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (equation < 0 || equation >= adjointer->nequations)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Invalid equation number %d.", equation);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (adjointer->callbacks.vec_destroy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_VEC_DESTROY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.vec_axpy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_VEC_AXPY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.mat_axpy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_MAT_AXPY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.mat_destroy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_MAT_DESTROY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }

  fwd_eqn = adjointer->equations[equation];
  *fwd_var = fwd_eqn.variable;
  ierr = adj_find_variable_data(&(adjointer->varhash), fwd_var, &fwd_data);
  assert(ierr == ADJ_OK);

  /* Check that we have all the forward values we need, before we start allocating stuff */
  for (i = 0; i < fwd_eqn.nblocks; i++)
  {
    adj_variable other_fwd_var;

    /* Check that we have the nonlinear block dependencies */
    if (fwd_eqn.blocks[i].has_nonlinear_block)
    {
      adj_nonlinear_block nl_block;

      nl_block = fwd_eqn.blocks[i].nonlinear_block;
      for (j=0; j < nl_block.ndepends; j++)
      {
        other_fwd_var = nl_block.depends[j];
        ierr = adj_has_variable_value(adjointer, other_fwd_var);
        if (ierr != ADJ_OK && !adj_variable_equal(fwd_var, &other_fwd_var, 1)) /* it's OK to not have the variable we're solving for, I suppose */
        {
          char buf[255];
          adj_variable_str(other_fwd_var, buf, 255);
          snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for variable %s, but don't have one.", buf);
          return adj_chkierr_auto(ADJ_ERR_NEED_VALUE);
        }
      }
    }

    /* Get the forward variable we want this to multiply */
    other_fwd_var = fwd_eqn.targets[i];
    if (adj_variable_equal(fwd_var, &other_fwd_var, 1)) continue; /* that term goes in the lhs */
    /* and now check it has a value */
    ierr = adj_has_variable_value(adjointer, other_fwd_var);
    if (ierr != ADJ_OK)
    {
      char buf[255];
      adj_variable_str(other_fwd_var, buf, 255);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for variable %s, but don't have one.", buf);
      return adj_chkierr_auto(ADJ_ERR_NEED_VALUE);
    }
  }

  /* Check the availability of the rhs dependencies */
  for (i=0; i < fwd_eqn.nrhsdeps; i++)
  {
    adj_variable other_fwd_var;

    other_fwd_var = fwd_eqn.rhsdeps[i];
    ierr = adj_has_variable_value(adjointer, other_fwd_var);
    if (ierr != ADJ_OK && !adj_variable_equal(fwd_var, &other_fwd_var, 1)) /* it's OK to not have the variable we're solving for, I suppose */
    {
      char buf[255];
      adj_variable_str(other_fwd_var, buf, 255);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for variable %s, but don't have one.", buf);
      return adj_chkierr_auto(ADJ_ERR_NEED_VALUE);
    }
  }

  /* --------------------------------------------------------------------------
   * Computation of A terms                                                  |
   * -------------------------------------------------------------------------- */

  /* fwd_data->targeting_equations what forward equations have nonzero blocks in the column of A associated with fwd_var. */

  {
    adj_block block;
    int blockcount = 0;
    for (i = 0; i < fwd_eqn.nblocks; i++)
    {
      if (adj_variable_equal(&(fwd_eqn.targets[i]), fwd_var, 1))
      {
        /* this is the right block */
        block = fwd_eqn.blocks[i];
        blockcount++;
        if (blockcount == 1) /* the first one we've found */
        {
          ierr = adj_evaluate_block_assembly(adjointer, block, lhs, rhs);
          if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        }
        else
        {
          adj_matrix lhs_tmp;
          adj_vector rhs_tmp;
          ierr = adj_evaluate_block_assembly(adjointer, block, &lhs_tmp, &rhs_tmp);
          if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
          adjointer->callbacks.vec_destroy(&rhs_tmp); /* we already have rhs from the first block assembly */
          adjointer->callbacks.mat_axpy(lhs, (adj_scalar) 1.0, lhs_tmp); /* add lhs_tmp to lhs */
          adjointer->callbacks.mat_destroy(&lhs_tmp);
        }
      }
    }
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
    assert(ierr == ADJ_OK); /* we should have them all, we checked for them earlier */

    ierr = adj_evaluate_block_action(adjointer, block, value, &rhs_tmp);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adjointer->callbacks.vec_axpy(rhs, (adj_scalar)-1.0, rhs_tmp);
    adjointer->callbacks.vec_destroy(&rhs_tmp);
  }

  /* And any forward source terms */
  if (fwd_eqn.rhs_callback != NULL) 
  {
    int has_output = -666;
    ierr = adj_evaluate_forward_source(adjointer, equation, &rhs_tmp, &has_output);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

    if (has_output)
    {
      adjointer->callbacks.vec_axpy(rhs, (adj_scalar)1.0, rhs_tmp);
      adjointer->callbacks.vec_destroy(&rhs_tmp);
    }
  }

  return ADJ_OK;

}

int adj_get_forward_solution(adj_adjointer* adjointer, int equation, adj_vector* soln, adj_variable* fwd_var)
{
  int ierr;
  adj_matrix lhs;
  adj_vector rhs;

  /* Check for the required callbacks */ 
  strncpy(adj_error_msg, "Need the solve data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.solve == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  ierr = adj_get_forward_equation(adjointer, equation, &lhs, &rhs, fwd_var);
  if (ierr != ADJ_OK)
    return adj_chkierr_auto(ierr);

  /* Solve the linear system */
  adjointer->callbacks.solve(*fwd_var, lhs, rhs, soln); 
  adjointer->callbacks.vec_destroy(&rhs);
  adjointer->callbacks.mat_destroy(&lhs);
  
  return ADJ_OK;
}

/* This function replays a forward solve starting from start_equation until (including) stop_equation */
/* To keep the memory foot print as small as possible, the routine forgets variables as soon as possible, */
/* keeping only variables which are needed to solve equations > stop_equation+1 */
/* If the checkpoint_last_timestep flag is ADJ_TRUE, then the solution variables of the equations with */
/* the same timestep as stop_equation will be  checkpointed. */
int adj_replay_forward_equations(adj_adjointer* adjointer, int start_equation, int stop_equation, int checkpoint_last_timestep)
{
  int equation, stop_timestep;
  int ierr;
  adj_vector soln;
  adj_variable var;
  adj_variable_data* var_data;
  adj_storage_data storage;


  /* Get the timstep of the last equation in the replay. Its solution will be recorded to memory */
  stop_timestep = adjointer->equations[stop_equation].variable.timestep;

  for (equation=start_equation; equation <= stop_equation; equation++)
  {
    /* We might have the solution of this equation already,
     * in which case we do not have to solve for it.
     */
    if ((adj_has_variable_value(adjointer, adjointer->equations[equation].variable) != ADJ_OK) ||
         (adjointer->revolve_data.overwrite == ADJ_TRUE))
    {
      if (adjointer->revolve_data.verbose)
        printf("Revolve: Replaying equation %i.\n", equation);
      ierr = adj_get_forward_solution(adjointer, equation, &soln, &var);
      if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

      /* Record the solution to memory. We always use adj_storage_memory_copy for recording as it is a safe choice */
      ierr = adj_storage_memory_copy(soln, &storage);
      if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      if (adjointer->revolve_data.overwrite)
      {
        ierr = adj_storage_set_overwrite(&storage, ADJ_TRUE);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        ierr = adj_storage_set_compare(&storage, ADJ_TRUE, adjointer->revolve_data.comparison_tolerance);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      }

      ierr = adj_record_variable(adjointer, var, storage);
      if (ierr<0)
        adj_chkierr(ierr);
      else if (ierr != ADJ_OK)
        return adj_chkierr_auto(ierr);

      adjointer->callbacks.vec_destroy(&soln);
    }
    else
    {
      var = adjointer->equations[equation].variable;
      if (adjointer->revolve_data.verbose)
        printf("Revolve: No need to replay equation %i.\n", equation);
    }

    /* Checkpoint the equation if desired */
    if (checkpoint_last_timestep == ADJ_TRUE && var.timestep == stop_timestep)
    {
      /* Mark the first equation of the last timestep as a checkpoint equation */
      if ((equation>0) && (adjointer->equations[equation-1].variable.timestep != adjointer->equations[equation].variable.timestep))
      {
        if (adjointer->revolve_data.verbose)
          printf("Revolve: Checkpoint equation %i in memory.\n", equation);
        ierr = adj_checkpoint_equation(adjointer, equation, ADJ_CHECKPOINT_STORAGE_MEMORY);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      }
      ierr = adj_find_variable_data(&(adjointer->varhash), &var, &var_data);
      assert(ierr == ADJ_OK);
      var_data->storage.storage_memory_is_checkpoint=ADJ_TRUE;
    }

    /* Forget everything that is not needed for future forward calculations */
    if (adjointer->revolve_data.overwrite != ADJ_TRUE)
    {
      ierr = adj_forget_forward_equation(adjointer, equation);
      if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    }
  }

  return ADJ_OK;
}

int adj_get_tlm_equation(adj_adjointer* adjointer, int equation, char* parameter, adj_matrix* lhs, adj_vector* rhs, adj_variable* tlm_var)
{
  int ierr;
  adj_equation fwd_eqn;
  adj_vector rhs_tmp;
  int i, j;
  adj_variable fwd_var;
  adj_variable_data* tlm_data;
  adj_variable_data* fwd_data;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING)
  {
    strncpy(adj_error_msg, "You have asked for a tangent linear model equation, but the adjointer has been deactivated.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (equation < 0 || equation >= adjointer->nequations)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Invalid equation number %d.", equation);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (adjointer->callbacks.vec_destroy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_VEC_DESTROY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.vec_axpy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_VEC_AXPY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.mat_axpy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_MAT_AXPY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.mat_destroy == NULL)
  {
    strncpy(adj_error_msg, "Need the ADJ_MAT_DESTROY_CB callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }

  fwd_eqn = adjointer->equations[equation];
  fwd_var = fwd_eqn.variable;
  memcpy(&fwd_var, &fwd_eqn.variable, sizeof(adj_variable));
  memcpy(tlm_var, &fwd_var, sizeof(adj_variable));
  tlm_var->type = ADJ_TLM; strncpy(tlm_var->functional, parameter, ADJ_NAME_LEN);

  /* Let's take care of the hash table. */
  ierr = adj_find_variable_data(&(adjointer->varhash), &fwd_var, &fwd_data);
  if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

  ierr = adj_find_variable_data(&(adjointer->varhash), tlm_var, &tlm_data);
  if (ierr == ADJ_ERR_HASH_FAILED)
  {
    /* It might not fail, if we have tried to fetch this equation already */

    ierr = adj_add_new_hash_entry(adjointer, tlm_var, &tlm_data);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

    tlm_data->equation = -2;
    tlm_data->ntargeting_equations = fwd_data->ntargeting_equations;
    tlm_data->targeting_equations = (int*) malloc(fwd_data->ntargeting_equations * sizeof(int));
    ADJ_CHKMALLOC(tlm_data->targeting_equations);
    memcpy(tlm_data->targeting_equations, fwd_data->targeting_equations, fwd_data->ntargeting_equations * sizeof(int));

    tlm_data->ndepending_equations = fwd_data->ndepending_equations;
    tlm_data->depending_equations = (int*) malloc(fwd_data->ndepending_equations * sizeof(int));
    ADJ_CHKMALLOC(tlm_data->depending_equations);
    memcpy(tlm_data->depending_equations, fwd_data->depending_equations, fwd_data->ndepending_equations * sizeof(int));

    tlm_data->nrhs_equations = fwd_data->nrhs_equations;
    tlm_data->rhs_equations = (int*) malloc(fwd_data->nrhs_equations * sizeof(int));
    ADJ_CHKMALLOC(tlm_data->rhs_equations);
    memcpy(tlm_data->rhs_equations, fwd_data->rhs_equations, fwd_data->nrhs_equations * sizeof(int));
  }


  /* Check that we have all the forward values we need, before we start allocating stuff */
  for (i = 0; i < fwd_eqn.nblocks; i++)
  {
    adj_variable other_fwd_var;

    /* Check that we have the nonlinear block dependencies */
    if (fwd_eqn.blocks[i].has_nonlinear_block)
    {
      adj_nonlinear_block nl_block;

      nl_block = fwd_eqn.blocks[i].nonlinear_block;
      for (j=0; j < nl_block.ndepends; j++)
      {
        other_fwd_var = nl_block.depends[j];
        ierr = adj_has_variable_value(adjointer, other_fwd_var);
        if (ierr != ADJ_OK && !adj_variable_equal(&fwd_var, &other_fwd_var, 1)) /* it's OK to not have the variable we're solving for, I suppose */
        {
          char buf[255];
          adj_variable_str(other_fwd_var, buf, 255);
          snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for variable %s, but don't have one.", buf);
          return adj_chkierr_auto(ADJ_ERR_NEED_VALUE);
        }
      }
    }

    /* Get the forward variable we want this to multiply */
    other_fwd_var = fwd_eqn.targets[i];
    if (adj_variable_equal(&fwd_var, &other_fwd_var, 1)) continue; /* that term goes in the lhs */
    other_fwd_var.type = ADJ_TLM;
    strncpy(other_fwd_var.functional, parameter, ADJ_NAME_LEN);
    /* and now check it has a value */
    ierr = adj_has_variable_value(adjointer, other_fwd_var);
    if (ierr != ADJ_OK)
    {
      char buf[255];
      adj_variable_str(other_fwd_var, buf, 255);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for variable %s, but don't have one.", buf);
      return adj_chkierr_auto(ADJ_ERR_NEED_VALUE);
    }
  }

  /* --------------------------------------------------------------------------
   * Computation of A terms                                                  |
   * -------------------------------------------------------------------------- */

  {
    adj_block block;
    int blockcount = 0;
    for (i = 0; i < fwd_eqn.nblocks; i++)
    {
      if (adj_variable_equal(&(fwd_eqn.targets[i]), &fwd_var, 1))
      {
        /* this is the right block */
        block = fwd_eqn.blocks[i];
        blockcount++;
        if (blockcount == 1) /* the first one we've found */
        {
          ierr = adj_evaluate_block_assembly(adjointer, block, lhs, rhs);
          if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        }
        else
        {
          adj_matrix lhs_tmp;
          adj_vector rhs_tmp;
          ierr = adj_evaluate_block_assembly(adjointer, block, &lhs_tmp, &rhs_tmp);
          if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
          adjointer->callbacks.vec_destroy(&rhs_tmp); /* we already have rhs from the first block assembly */
          adjointer->callbacks.mat_axpy(lhs, (adj_scalar) 1.0, lhs_tmp); /* add lhs_tmp to lhs */
          adjointer->callbacks.mat_destroy(&lhs_tmp);
        }
      }
    }
  }

  /* Great! Now let's assemble the RHS contributions of A. */

  /* Now loop through the off-diagonal blocks of A. */
  for (i = 0; i < fwd_eqn.nblocks; i++)
  {
    adj_block block;
    adj_variable tlm_var;
    adj_vector value;

    /* Ignore the diagonal block */
    if (adj_variable_equal(&fwd_eqn.targets[i], &fwd_var, 1))
      continue;

    /* Get the TLM variable we want this block to multiply with */
    tlm_var = fwd_eqn.targets[i];
    tlm_var.type = ADJ_TLM;
    strncpy(tlm_var.functional, parameter, ADJ_NAME_LEN);

    block = fwd_eqn.blocks[i];

    /* and now get its value */
    ierr = adj_get_variable_value(adjointer, tlm_var, &value);
    assert(ierr == ADJ_OK); /* we should have them all, we checked for them earlier */

    ierr = adj_evaluate_block_action(adjointer, block, value, &rhs_tmp);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    adjointer->callbacks.vec_axpy(rhs, (adj_scalar)-1.0, rhs_tmp);
    adjointer->callbacks.vec_destroy(&rhs_tmp);
  }

  /* --------------------------------------------------------------------------
   * Computation of G terms                                                  |
   * -------------------------------------------------------------------------- */

  /* We need to loop through the dependencies of this equation; each one of those
     dependencies will produce a G entry. */

  {
    int nderivs; /* these two are the raw derivatives to compute */
    adj_nonlinear_block_derivative* derivs;

    int nnew_derivs; /* and these two are after derivative simplification */
    adj_nonlinear_block_derivative* new_derivs;
    int l, k;

    /* First, we find how many derivatives we're going to need, so that we can malloc
       appropriately. */
    nderivs = 0;
    for (i = 0; i < fwd_eqn.nblocks; i++)
    {
      if (fwd_eqn.blocks[i].has_nonlinear_block)
      {
        nderivs += fwd_eqn.blocks[i].nonlinear_block.ndepends;
      }
    }

    /* Now that we have nderivs, let's use it */

    if (nderivs > 0)
    {
      derivs = (adj_nonlinear_block_derivative*) malloc(nderivs * sizeof(adj_nonlinear_block_derivative));
      ADJ_CHKMALLOC(derivs);
      l = 0;
      for (i = 0; i < fwd_eqn.nblocks; i++)
      {
        if (fwd_eqn.blocks[i].has_nonlinear_block)
        {
          for (k = 0; k < fwd_eqn.blocks[i].nonlinear_block.ndepends; k++)
          {
            adj_vector target;
            ierr = adj_get_variable_value(adjointer, fwd_eqn.targets[i], &target);
            if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

            ierr = adj_create_nonlinear_block_derivative(adjointer, fwd_eqn.blocks[i].nonlinear_block, fwd_eqn.blocks[i].coefficient, fwd_eqn.blocks[i].nonlinear_block.depends[k], target, fwd_eqn.blocks[i].hermitian, &derivs[l]);
            l++;
          }
        }
      }

      /* OK, Here's where we do our simplifications; this can be a significant optimisation */
      /* .......................................................................................... */
      ierr = adj_simplify_derivatives(adjointer, nderivs, derivs, &nnew_derivs, &new_derivs);
      for (l = 0; l < nderivs; l++)
      {
        ierr = adj_destroy_nonlinear_block_derivative(adjointer, &derivs[l]);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      }
      free(derivs);

      /* Now, we go and evaluate each one of the derivatives, assembling or acting as appropriate */
      for (l = 0; l < nnew_derivs; l++)
      {
        if (adj_variable_equal(&new_derivs[l].variable, &fwd_var, 1))
        {
          /* This G-block is on the diagonal, and so we must assemble it .. later */
          snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Sorry, we can't handle G-blocks on the diagonal (yet).");
          return adj_chkierr_auto(ADJ_ERR_NOT_IMPLEMENTED);
        }
        else
        {
          /* This G-block is NOT on the diagonal, so we only need its action */
          adj_variable tlm_associated;
          adj_vector tlm_value;

          tlm_associated = new_derivs[l].variable;
          tlm_associated.type = ADJ_TLM;
          strncpy(tlm_associated.functional, parameter, ADJ_NAME_LEN);
          ierr = adj_get_variable_value(adjointer, tlm_associated, &tlm_value);
          if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

          /* And now we are ready */
          ierr = adj_evaluate_nonlinear_derivative_action(adjointer, 1, &new_derivs[l], tlm_value, rhs);
          if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        }
      }

      for (l = 0; l < nnew_derivs; l++)
      {
        ierr = adj_destroy_nonlinear_block_derivative(adjointer, &new_derivs[l]);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
      }
      free(new_derivs);
    }
  }

  /* --------------------------------------------------------------------------
   * Computation of R terms                                                  |
   * -------------------------------------------------------------------------- */

  {
    for (i = 0; i < fwd_eqn.nrhsdeps; i++)
    {
      /* Does this R contribute to the adjoint matrix ... */
      if (adj_variable_equal(&fwd_var, &fwd_eqn.rhsdeps[i], 1))
      {
        adj_matrix r;
        ierr = adj_evaluate_rhs_deriv_assembly(adjointer, fwd_eqn, ADJ_FALSE, &r);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
        adjointer->callbacks.mat_axpy(lhs, (adj_scalar) -1.0, r); /* Subtract the R contribution from the adjoint lhs */
        adjointer->callbacks.mat_destroy(&r);
      }
      /* ... or to the right-hand side of the adjoint system? */
      else
      {
        adj_vector deriv_action;
        int has_output;

        has_output = -666;

        adj_variable contraction_var;
        adj_vector contraction;

        contraction_var = fwd_eqn.rhsdeps[i];
        contraction_var.type = ADJ_TLM;
        strncpy(contraction_var.functional, parameter, ADJ_NAME_LEN);
        ierr = adj_get_variable_value(adjointer, contraction_var, &contraction);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        ierr = adj_evaluate_rhs_deriv_action(adjointer, fwd_eqn, fwd_eqn.rhsdeps[i], contraction, ADJ_FALSE, &deriv_action, &has_output);
        if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);

        if (has_output)
        {
          /* Now that we have the contribution, we need to add it to the adjoint right hand side */
          adjointer->callbacks.vec_axpy(rhs, (adj_scalar)1.0, deriv_action);
          adjointer->callbacks.vec_destroy(&deriv_action);
        }
      }
    }
  }

  /* And any tangent linear source terms */
  {
    adj_vector rhs_tmp;
    int has_psrc;
    has_psrc = -666;
    ierr = adj_evaluate_parameter_source(adjointer, fwd_var, parameter, &rhs_tmp, &has_psrc);
    if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    if (has_psrc)
    {
      adjointer->callbacks.vec_axpy(rhs, (adj_scalar)1.0, rhs_tmp);
      adjointer->callbacks.vec_destroy(&rhs_tmp);
    }
  }

  return ADJ_OK;

}

int adj_get_tlm_solution(adj_adjointer* adjointer, int equation, char* parameter, adj_vector* soln, adj_variable* tlm_var)
{
  int ierr;
  adj_matrix lhs;
  adj_vector rhs;

  /* Check for the required callbacks */ 
  strncpy(adj_error_msg, "Need the solve data callback, but it hasn't been supplied.", ADJ_ERROR_MSG_BUF);
  if (adjointer->callbacks.solve == NULL) return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  strncpy(adj_error_msg, "", ADJ_ERROR_MSG_BUF);

  ierr = adj_get_tlm_equation(adjointer, equation, parameter, &lhs, &rhs, tlm_var);
  if (ierr != ADJ_OK)
    return adj_chkierr_auto(ierr);

  /* Solve the linear system */
  adjointer->callbacks.solve(*tlm_var, lhs, rhs, soln); 
  adjointer->callbacks.vec_destroy(&rhs);
  adjointer->callbacks.mat_destroy(&lhs);

  return ADJ_OK;
}
