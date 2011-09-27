#include "libadjoint/adj_adjointer_routines.h"

int adj_create_adjointer(adj_adjointer* adjointer)
{
  int i;

  adjointer->nequations = 0;
  adjointer->equations_sz = 0;
  adjointer->equations = NULL;

  adjointer->ntimesteps = 0;
  adjointer->timestep_data = NULL;

  adjointer->varhash = NULL;
  adjointer->vardata.firstnode = NULL;
  adjointer->vardata.lastnode = NULL;

  adjointer->callbacks.vec_duplicate = NULL;
  adjointer->callbacks.vec_axpy = NULL;
  adjointer->callbacks.vec_destroy = NULL;
  adjointer->callbacks.vec_set_values = NULL;
  adjointer->callbacks.vec_get_size = NULL;
  adjointer->callbacks.vec_divide = NULL;
  adjointer->callbacks.vec_get_norm = NULL;
  adjointer->callbacks.vec_set_random = NULL;
  adjointer->callbacks.vec_dot_product = NULL;
  adjointer->callbacks.vec_to_file = NULL;
  adjointer->callbacks.vec_from_file = NULL;

  adjointer->callbacks.mat_duplicate = NULL;
  adjointer->callbacks.mat_axpy = NULL;
  adjointer->callbacks.mat_destroy = NULL;

  adjointer->callbacks.solve = NULL;

  adjointer->revolve_data.steps = 0;
  adjointer->revolve_data.snaps = 0;
  adjointer->revolve_data.snaps_in_ram = 0;
  adjointer->revolve_data.revolve.ptr = NULL;
  adjointer->revolve_data.current_action = -1;
  adjointer->revolve_data.current_timestep = -999999999;

  adjointer->nonlinear_colouring_list.firstnode = NULL;
  adjointer->nonlinear_colouring_list.lastnode = NULL;
  adjointer->nonlinear_action_list.firstnode = NULL;
  adjointer->nonlinear_action_list.lastnode = NULL;
  adjointer->nonlinear_derivative_action_list.firstnode = NULL;
  adjointer->nonlinear_derivative_action_list.lastnode = NULL;
  adjointer->nonlinear_derivative_assembly_list.firstnode = NULL;
  adjointer->nonlinear_derivative_assembly_list.lastnode = NULL;
  adjointer->block_action_list.firstnode = NULL;
  adjointer->block_action_list.lastnode = NULL;
  adjointer->block_assembly_list.firstnode = NULL;
  adjointer->block_assembly_list.lastnode = NULL;
  adjointer->functional_list.firstnode = NULL;
  adjointer->functional_list.lastnode = NULL;
  adjointer->functional_derivative_list.firstnode = NULL;
  adjointer->functional_derivative_list.lastnode = NULL;
  adjointer->forward_source_callback = NULL;

  for (i = 0; i < ADJ_NO_OPTIONS; i++)
    adjointer->options[i] = 0; /* 0 is the default for all options */

  return ADJ_OK;
}

int adj_destroy_adjointer(adj_adjointer* adjointer)
{
  int i;
  int ierr;
  adj_variable_data* data_ptr;
  adj_variable_data* data_ptr_tmp;
  adj_op_callback* cb_ptr;
  adj_op_callback* cb_ptr_tmp;
  adj_func_callback* func_cb_ptr;
  adj_func_callback* func_cb_ptr_tmp;
  adj_func_deriv_callback* func_deriv_cb_ptr;
  adj_func_deriv_callback* func_deriv_cb_ptr_tmp;
  adj_functional_data* functional_data_ptr_next = NULL;
  adj_functional_data* functional_data_ptr = NULL;

  for (i = 0; i < adjointer->nequations; i++)
  {
    ierr = adj_destroy_equation(&(adjointer->equations[i]));
    if (ierr != ADJ_OK) return ierr;
  }
  if (adjointer->equations != NULL) free(adjointer->equations);

  if (adjointer->timestep_data != NULL)
  {
    for (i = 0; i < adjointer->ntimesteps; i++)
      functional_data_ptr = adjointer->timestep_data[i].functional_data_start;
      while (functional_data_ptr != NULL)
      {
        functional_data_ptr_next = functional_data_ptr->next;
        if (functional_data_ptr->dependencies != NULL) free(functional_data_ptr->dependencies);
        free(functional_data_ptr);
        functional_data_ptr = functional_data_ptr_next;
      }
    free(adjointer->timestep_data);
  }

  data_ptr = adjointer->vardata.firstnode;
  while (data_ptr != NULL)
  {
    ierr = adj_destroy_variable_data(adjointer, data_ptr);
    if (ierr != ADJ_OK) return ierr;
    data_ptr_tmp = data_ptr;
    data_ptr = data_ptr->next;
    free(data_ptr_tmp);
  }

  cb_ptr = adjointer->nonlinear_colouring_list.firstnode;
  while(cb_ptr != NULL)
  {
    cb_ptr_tmp = cb_ptr;
    cb_ptr = cb_ptr->next;
    free(cb_ptr_tmp);
  }

  cb_ptr = adjointer->nonlinear_action_list.firstnode;
  while(cb_ptr != NULL)
  {
    cb_ptr_tmp = cb_ptr;
    cb_ptr = cb_ptr->next;
    free(cb_ptr_tmp);
  }

  cb_ptr = adjointer->nonlinear_derivative_action_list.firstnode;
  while(cb_ptr != NULL)
  {
    cb_ptr_tmp = cb_ptr;
    cb_ptr = cb_ptr->next;
    free(cb_ptr_tmp);
  }

  cb_ptr = adjointer->nonlinear_derivative_assembly_list.firstnode;
  while(cb_ptr != NULL)
  {
    cb_ptr_tmp = cb_ptr;
    cb_ptr = cb_ptr->next;
    free(cb_ptr_tmp);
  }

  cb_ptr = adjointer->block_action_list.firstnode;
  while(cb_ptr != NULL)
  {
    cb_ptr_tmp = cb_ptr;
    cb_ptr = cb_ptr->next;
    free(cb_ptr_tmp);
  }

  cb_ptr = adjointer->block_assembly_list.firstnode;
  while(cb_ptr != NULL)
  {
    cb_ptr_tmp = cb_ptr;
    cb_ptr = cb_ptr->next;
    free(cb_ptr_tmp);
  }

  func_cb_ptr = adjointer->functional_list.firstnode;
  while(func_cb_ptr != NULL)
  {
    func_cb_ptr_tmp = func_cb_ptr;
    func_cb_ptr = func_cb_ptr->next;
    free(func_cb_ptr_tmp);
  }
  
  func_deriv_cb_ptr = adjointer->functional_derivative_list.firstnode;
  while(func_deriv_cb_ptr != NULL)
  {
    func_deriv_cb_ptr_tmp = func_deriv_cb_ptr;
    func_deriv_cb_ptr = func_deriv_cb_ptr->next;
    free(func_deriv_cb_ptr_tmp);
  }

  adj_create_adjointer(adjointer);
  return ADJ_OK;
}

int adj_deactivate_adjointer(adj_adjointer* adjointer)
{
  return adj_set_option(adjointer, ADJ_ACTIVITY, ADJ_ACTIVITY_NOTHING);
}

int adj_set_checkpoint_strategy(adj_adjointer* adjointer, int strategy)
{
  return adj_set_option(adjointer, ADJ_CHECKPOINT_STRATEGY, strategy);
}

int adj_get_checkpoint_strategy(adj_adjointer* adjointer, int* strategy)
{
  *strategy = adjointer->options[ADJ_CHECKPOINT_STRATEGY];
  return ADJ_OK;
}

int adj_set_revolve_options(adj_adjointer* adjointer, int steps, int snaps, int snaps_in_ram)
{
  adjointer->revolve_data.steps=steps;  
  adjointer->revolve_data.snaps=snaps;  
  adjointer->revolve_data.snaps_in_ram=snaps_in_ram;  
  return ADJ_OK;
}

int adj_register_equation(adj_adjointer* adjointer, adj_equation equation, int* checkpoint_storage)
{
  adj_variable_data* data_ptr;
  int ierr;
  int i;
  int j;
  int cs; /* checkpoint strategy */

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;

  /* Let's check we haven't solved for this variable before */
  ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.variable), &data_ptr);
  if (ierr != ADJ_ERR_HASH_FAILED)
  {
    /* We may have legitimately seen this before, if it's been registered as a functional dependency.
       So we have to check its hash table entry for the equation */
    if (data_ptr->equation >= 0)
    {
      char buf[ADJ_NAME_LEN];
      adj_variable_str(equation.variable, buf, ADJ_NAME_LEN);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "We have already registered an equation for variable %s.", buf);
      return ADJ_ERR_INVALID_INPUTS;
    }
  }

  /* Let's check the timesteps match up */
  if (adjointer->nequations == 0) /* we haven't registered any equations yet */
  {
    if (equation.variable.timestep != 0) /* this isn't timestep 0 */
    {
      strncpy(adj_error_msg, "The first equation registered must have timestep 0.", ADJ_ERROR_MSG_BUF);
      return ADJ_ERR_INVALID_INPUTS;
    }
  }
  else /* we have registered an equation before */
  {
    int old_timestep;
    /* if  (not same timestep as before)                     &&  (not the next timestep) */
    old_timestep = adjointer->equations[adjointer->nequations-1].variable.timestep;
    if ((equation.variable.timestep != old_timestep) && (equation.variable.timestep != old_timestep+1))
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, \
          "Timestep numbers must either stay the same or increment by one. Valid values are %d or %d, but you have supplied %d.", \
          old_timestep, old_timestep+1, equation.variable.timestep);
      return ADJ_ERR_INVALID_INPUTS;
    }
  }

  if (equation.variable.auxiliary)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(equation.variable, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Cannot register an equation for an auxiliary variable %s.", buf);
  }

  /* check that any targets that aren't auxiliary and aren't the subject of the current equation
     have equations registered for them (i.e., the system is lower-triangular) */
  for (i = 0; i < equation.nblocks; i++)
  {
    if (!equation.targets[i].auxiliary && !adj_variable_equal(&equation.targets[i], &equation.variable, 1))
    {
      adj_variable_data* hash_ptr;
      int ierr;
      ierr = adj_find_variable_data(&adjointer->varhash, &equation.targets[i], &hash_ptr);
      if (ierr == ADJ_ERR_HASH_FAILED)
      {
        char buf[ADJ_NAME_LEN];
        adj_variable_str(equation.targets[i], buf, ADJ_NAME_LEN);
        snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Have %s as a target in an equation, but have never seen that variable before.", buf);
        return ierr;
      }

      if (hash_ptr->equation < 0)
      {
        char buf[ADJ_NAME_LEN];
        adj_variable_str(equation.targets[i], buf, ADJ_NAME_LEN);
        snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Have %s as a target in an equation, but that variable has not been solved for previously.", buf);
        return ADJ_ERR_INVALID_INPUTS;
      }
    }
  }

  /* OK. We're good to go. */

  /* Let's add it to the hash table. */
  /* ierr is from the call to adj_find_variable_data. If it is ADJ_ERR_HASH_FAILED, it means
     we have to add it in to the hash table. */
  if (ierr == ADJ_ERR_HASH_FAILED)
  {
    ierr = adj_add_new_hash_entry(adjointer, &(equation.variable), &data_ptr);
    if (ierr != ADJ_OK) return ierr;
  }
  data_ptr->equation = adjointer->nequations;
  /* OK. Next create an entry for the adj_equation in the adjointer. */

  /* Check we have enough room, and if not, make some */
  if (adjointer->nequations == adjointer->equations_sz)
  {
    adjointer->equations = (adj_equation*) realloc(adjointer->equations, (adjointer->equations_sz + ADJ_PREALLOC_SIZE) * sizeof(adj_equation));
    ADJ_CHKMALLOC(adjointer->equations);
    adjointer->equations_sz = adjointer->equations_sz + ADJ_PREALLOC_SIZE;
  }

  adjointer->nequations++;
  adjointer->equations[adjointer->nequations - 1] = equation;

  /* Do any necessary recording of timestep indices */
  if (adjointer->ntimesteps < equation.variable.timestep + 1) /* adjointer->ntimesteps should be at least equation.variable.timestep + 1 */
  {
    ierr = adj_extend_timestep_data(adjointer, equation.variable.timestep + 1); /* extend the array as necessary */
    if (ierr != ADJ_OK) return ierr;
  }
  if (adjointer->timestep_data[equation.variable.timestep].start_equation == -1) /* -1 is the sentinel value for unset */
  {
    adjointer->timestep_data[equation.variable.timestep].start_equation = adjointer->nequations - 1; /* fill in the start equation */
  }

  /* now we have copies of the pointer to the arrays of targets, blocks, rhs deps. */
  /* but for consistency, any libadjoint object that the user creates, he must destroy --
     it's simpler that way. */
  /* so we're going to make our own copies, so that the user can destroy his. */

  /* blocks */
  adjointer->equations[adjointer->nequations - 1].blocks = (adj_block*) malloc(equation.nblocks * sizeof(adj_block));
  ADJ_CHKMALLOC(adjointer->equations[adjointer->nequations - 1].blocks);
  memcpy(adjointer->equations[adjointer->nequations - 1].blocks, equation.blocks, equation.nblocks * sizeof(adj_block));
  for (i = 0; i < equation.nblocks; i++)
  {
    if (equation.blocks[i].has_nonlinear_block)
    {
      int ierr;
      ierr = adj_copy_nonlinear_block(equation.blocks[i].nonlinear_block, &adjointer->equations[adjointer->nequations - 1].blocks[i].nonlinear_block);
      if (ierr != ADJ_OK) return ierr;
    }
  }

  /* targets */
  adjointer->equations[adjointer->nequations - 1].targets = (adj_variable*) malloc(equation.nblocks * sizeof(adj_variable));
  ADJ_CHKMALLOC(adjointer->equations[adjointer->nequations - 1].targets);
  memcpy(adjointer->equations[adjointer->nequations - 1].targets, equation.targets, equation.nblocks * sizeof(adj_variable));
  if (equation.nrhsdeps > 0)
  {
    adjointer->equations[adjointer->nequations - 1].rhsdeps = (adj_variable*) malloc(equation.nrhsdeps * sizeof(adj_variable));
    ADJ_CHKMALLOC(adjointer->equations[adjointer->nequations - 1].rhsdeps);
    memcpy(adjointer->equations[adjointer->nequations - 1].rhsdeps, equation.rhsdeps, equation.nrhsdeps * sizeof(adj_variable));
  }

  /* Now find all the entries we need to update in the hash table, and update them */
  /* First: targeting equations */

  for (i = 0; i < equation.nblocks; i++)
  {
    ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.targets[i]), &data_ptr);
    if (ierr != ADJ_OK)
    {
      char buf[ADJ_NAME_LEN];
      adj_variable_str(equation.targets[i], buf, ADJ_NAME_LEN);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The equation to be registered has a block is targeting %s, but I do not have an equation for that variable yet.", buf);
      return ierr;
    }

    /* this is already guaranteed to be a unique entry -- we have never seen this equation before.
       so we don't need adj_append_unique */
    data_ptr->ntargeting_equations++;
    data_ptr->targeting_equations = (int*) realloc(data_ptr->targeting_equations, data_ptr->ntargeting_equations * sizeof(int));
    ADJ_CHKMALLOC(data_ptr->targeting_equations);
    data_ptr->targeting_equations[data_ptr->ntargeting_equations - 1] = adjointer->nequations - 1;
  }

  /* Next: nonlinear dependencies of the right hand side */
  for (i = 0; i < equation.nrhsdeps; i++)
  {
    ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.rhsdeps[i]), &data_ptr);
    if (ierr == ADJ_ERR_HASH_FAILED && equation.rhsdeps[i].auxiliary)
    {
      /* It's ok if it's auxiliary -- it legitimately can be the first time we've seen it */
      adj_variable_data* new_data;
      ierr = adj_add_new_hash_entry(adjointer, &(equation.rhsdeps[i]), &new_data);
      if (ierr != ADJ_OK) return ierr;
      new_data->equation = -1; /* it doesn't have an equation */
      data_ptr = new_data;
    }
    else
    {
      return ierr;
    }

    /* Now data_ptr points to the data we're storing */
    ierr = adj_append_unique(&(data_ptr->rhs_equations), &(data_ptr->nrhs_equations), adjointer->nequations - 1);
    if (ierr != ADJ_OK) return ierr;
  }

  /* Now we need to record what we need for which adjoint equation from the rhs dependencies */
  /* We do this after the previous loop to make sure that all of the dependencies exist in the hash table. */
  for (i = 0; i < equation.nrhsdeps; i++)
  {
    /* We will loop through all other dependencies and register that the adjoint equation
       associated with equation.rhsdeps[i] depends on all other equation.rhsdeps[j] */
    int eqn_no;

    if (equation.rhsdeps[i].auxiliary)
    {
      /* We don't solve any equation for it, and auxiliary variables thus don't have any associated
         adjoint equation -- so we don't need to register any dependencies for it */
      continue;
    }

    ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.rhsdeps[i]), &data_ptr);
    if (ierr != ADJ_OK) return ierr;

    eqn_no = data_ptr->equation;
    assert(eqn_no >= 0);

    for (j = 0; j < equation.nrhsdeps; j++)
    {
      ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.rhsdeps[j]), &data_ptr);
      if (ierr != ADJ_OK) return ierr;
      ierr = adj_append_unique(&(data_ptr->adjoint_equations), &(data_ptr->nadjoint_equations), eqn_no); /* dependency j is necessary for equation i */
      if (ierr != ADJ_OK) return ierr;
    }
  }

  /* And finally nonlinear dependencies of the left hand side */

  for (i = 0; i< equation.nblocks; i++)
  {
    if (equation.blocks[i].has_nonlinear_block)
    {
      for (j = 0; j < equation.blocks[i].nonlinear_block.ndepends; j++)
      {
        /* Register that this equation depends on this variable */
        ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.blocks[i].nonlinear_block.depends[j]), &data_ptr);
        if (ierr == ADJ_ERR_HASH_FAILED && equation.blocks[i].nonlinear_block.depends[j].auxiliary)
        {
          /* It's ok if it's auxiliary -- it legitimately can be the first time we've seen it */
          adj_variable_data* new_data;
          ierr = adj_add_new_hash_entry(adjointer, &(equation.blocks[i].nonlinear_block.depends[j]), &new_data);
          if (ierr != ADJ_OK) return ierr;
          new_data->equation = -1; /* it doesn't have an equation */
          data_ptr = new_data;
        }
        else if (ierr == ADJ_ERR_HASH_FAILED)
        {
          return ierr;
        }
        ierr = adj_append_unique(&(data_ptr->depending_equations), &(data_ptr->ndepending_equations), adjointer->nequations - 1);
        if (ierr != ADJ_OK) return ierr;
      }
    }
  }

  /* And now perform the updates for .adjoint_equations implied by the existence of these dependencies. 
     This is probably the hardest, most mindbending thing in the whole library -- sorry. There's no way
     to really make this easy; you just have to work through it. */
  for (i = 0; i < equation.nblocks; i++)
  {
    if (equation.blocks[i].has_nonlinear_block)
    {
      adj_variable_data* block_target_data; /* fetch the hash entry associated with the target of this block */
      ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.targets[i]), &block_target_data);
      if (ierr != ADJ_OK) return ierr;

      for (j = 0; j < equation.blocks[i].nonlinear_block.ndepends; j++)
      {
        int k;
        adj_variable_data* j_data;

        /* j_data ALWAYS refers to the data associated with the j'th dependency, throughout this whole loop */
        ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.blocks[i].nonlinear_block.depends[j]), &j_data);
        if (ierr != ADJ_OK) return ierr;

        /* One set of dependencies: the (adjoint equation of) (the target of this block) (needs) (this dependency) */
        ierr = adj_append_unique(&(j_data->adjoint_equations), &(j_data->nadjoint_equations), block_target_data->equation);
        if (ierr != ADJ_OK) return ierr;

        /* Another set of dependencies: the (adjoint equation of) (the j'th dependency) (needs) (the target of this block) */
        ierr = adj_append_unique(&(block_target_data->adjoint_equations), &(block_target_data->nadjoint_equations), j_data->equation);
        if (ierr != ADJ_OK) return ierr;

        /* Now we loop over all the dependencies again and fill in the cross-dependencies */
        for (k = 0; k < equation.blocks[i].nonlinear_block.ndepends; k++)
        {
          adj_variable_data* k_data;

          /* k_data ALWAYS refers to the data associated with the k'th dependency, throughout this whole loop */
          ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.blocks[i].nonlinear_block.depends[k]), &k_data);
          if (ierr != ADJ_OK) return ierr;

          /* Another set of dependencies: the (adjoint equation of) (the j'th dependency) (needs) (the k'th dependency) */
          ierr = adj_append_unique(&(k_data->adjoint_equations), &(k_data->nadjoint_equations), j_data->equation);
          if (ierr != ADJ_OK) return ierr;
        }
      }
    }
  }

  /* Set the checkpoint flag */
  /* TODO: Do all this logic in adj_create_equation and set the flag in adj_create_equation??? */
  *checkpoint_storage = ADJ_CHECKPOINT_STORAGE_NONE;

  ierr = adj_get_checkpoint_strategy(adjointer, &cs);
  if (ierr != ADJ_OK) return ierr;

  if ((cs==ADJ_CHECKPOINT_REVOLVE_OFFLINE) || (cs==ADJ_CHECKPOINT_REVOLVE_MULTISTAGE) || (cs==ADJ_CHECKPOINT_REVOLVE_ONLINE))
  {
    ierr = adj_get_revolve_checkpoint_storage(adjointer, equation, checkpoint_storage);
    if (ierr != ADJ_OK) return ierr;
    adjointer->equations[adjointer->nequations - 1].checkpoint_type = *checkpoint_storage;
  }

  return ADJ_OK;
}

int adj_get_revolve_checkpoint_storage(adj_adjointer* adjointer, adj_equation equation, int *checkpoint_storage)
{
  int cs, ierr;
  char buf[ADJ_NAME_LEN];
  int oldcapo, capo;

  ierr = adj_get_checkpoint_strategy(adjointer, &cs);
  if (ierr != ADJ_OK) return ierr;

  /* Initialise Revolve if not done before and otherwise do some consistency checks*/
  if (adjointer->revolve_data.revolve.ptr != NULL)
  {
    if ((adjointer->revolve_data.current_action != CACTION_ADVANCE) && (adjointer->revolve_data.current_action != CACTION_FIRSTRUN))
    {
      char buf[ADJ_NAME_LEN];
      adj_variable_str(equation.variable, buf, ADJ_NAME_LEN);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The adjointer and revolve are not consistent (in adj_register_equation of variable %s). This can happen when one tries to register an equation after solving the first adjoint equation.", buf);
      return ADJ_ERR_REVOLVE_ERROR;
    }
  }
  else
  {
    ierr = adj_initialise_revolve(adjointer);
    if (ierr != ADJ_OK) return ierr;

    /* Set the intial revolve state */
    adjointer->revolve_data.current_action = revolve(adjointer->revolve_data.revolve);
    adjointer->revolve_data.current_timestep = equation.variable.timestep;
    if (equation.variable.timestep!=0) 
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "With revolve as checkpoint strategy the first equation has solve for a variable at timestep 0.");
      return ADJ_ERR_INVALID_INPUTS;
    }
  }

  /* Check that the equations are registered chronologically */
  if ((adjointer->revolve_data.current_timestep != equation.variable.timestep) && (adjointer->revolve_data.current_timestep+1 != equation.variable.timestep))
  {
    adj_variable_str(equation.variable, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "With revolve as checkpoint strategy the equations have to be registered chronologically, but equation for variable %s is out of order.", buf);
    return ADJ_ERR_INVALID_INPUTS;
  }

  /* Update revolve state */
  adjointer->revolve_data.current_timestep = equation.variable.timestep;

  /* Determine if a checkpoint is requested and of which type */
  switch (adjointer->revolve_data.current_action)
  {

    case CACTION_ADVANCE:
      capo = revolve_getcapo(adjointer->revolve_data.revolve);
      oldcapo = revolve_getoldcapo(adjointer->revolve_data.revolve);

      /* make sure that Revolve and the adjointer are in sync */
      if ((adjointer->revolve_data.current_timestep <= oldcapo) || (adjointer->revolve_data.current_timestep > capo))
      {
        adj_variable_str(equation.variable, buf, ADJ_NAME_LEN);
        snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "An internal error occured: The adjointer and revolve are out of sync (in adj_register_equation of variable %s).", buf);
        return ADJ_ERR_REVOLVE_ERROR;
      }

      if (adjointer->revolve_data.current_timestep == capo)
        adjointer->revolve_data.current_action = revolve(adjointer->revolve_data.revolve);

      /* In the case we want to take a checkpoint, we do not break here and execute the next case as well */
      if (adjointer->revolve_data.current_action != CACTION_TAKESHOT)
        break;

      case CACTION_TAKESHOT:
        if (cs==ADJ_CHECKPOINT_REVOLVE_MULTISTAGE)
        {
          if (revolve_getwhere(adjointer->revolve_data.revolve))
              *checkpoint_storage = ADJ_CHECKPOINT_STORAGE_MEMORY;
          else
              *checkpoint_storage = ADJ_CHECKPOINT_STORAGE_DISK;
        }
        else
          *checkpoint_storage = ADJ_CHECKPOINT_STORAGE_DISK;

        adjointer->revolve_data.current_action = revolve(adjointer->revolve_data.revolve);
        break;

      case CACTION_FIRSTRUN:
        break;

      case CACTION_ERROR:
        adj_variable_str(equation.variable, buf, ADJ_NAME_LEN);
        snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "An internal error occured: Irregular termination of revolve (in adj_register_equation of variable %s).", buf);
        return ADJ_ERR_REVOLVE_ERROR;

      default: /* This case includes CACTION_YOUTURN which is only expected when restoring from a checkpoint */
        adj_variable_str(equation.variable, buf, ADJ_NAME_LEN);
        snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "An internal error occured: The adjointer and revolve are out of sync (in adj_register_equation of variable %s).", buf);
        return ADJ_ERR_REVOLVE_ERROR;
  }

  return ADJ_OK;
}

int adj_initialise_revolve(adj_adjointer* adjointer)
{
  int steps = adjointer->revolve_data.steps;
  int snaps = adjointer->revolve_data.snaps;
  int snaps_in_ram = adjointer->revolve_data.snaps_in_ram;
  int cs, ierr;

  ierr = adj_get_checkpoint_strategy(adjointer, &cs);
  if (ierr != ADJ_OK) return ierr;

  /* Offline checkpointing */
  if (cs==ADJ_CHECKPOINT_REVOLVE_OFFLINE) 
  {
    if ((steps>0) && (snaps>0) && (snaps_in_ram<=0))
      adjointer->revolve_data.revolve = revolve_create_offline(steps, snaps);
    else
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You chose to use offline revolve as checkpointing strategy but have not configured it correctly. Make sure you call adj_set_revolve_options with positive 'steps' and 'snaps' arguments.");
      return ADJ_ERR_INVALID_INPUTS;
    }
  }

  /* Offline checkpointing with different stores */
  else if (cs==ADJ_CHECKPOINT_REVOLVE_MULTISTAGE)
  { 
    if ((steps>0) && (snaps>0) && (snaps_in_ram>0))
      adjointer->revolve_data.revolve = revolve_create_multistage(steps, snaps, snaps_in_ram);
    else
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You chose to use multistage revolve as checkpointing strategy but have not configured it correctly. Make sure you call adj_set_revolve_options with positive 'steps', 'snaps' and 'snaps_in_ram' arguments.");
      return ADJ_ERR_INVALID_INPUTS;
    }
  }

  /* Online checkpointing */
  else if (cs==ADJ_CHECKPOINT_REVOLVE_ONLINE) 
  {
    if (snaps>0)
      adjointer->revolve_data.revolve = revolve_create_online(snaps);
    else
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You chose to use online revolve as checkpointing strategy but have not configured it correctly. Make sure you call adj_set_revolve_options with a positive 'snaps' argument.");
      return ADJ_ERR_INVALID_INPUTS;
    }
  }

  else
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You chose to use revolve as checkpointing strategy but have not configured it correctly.");
    return ADJ_ERR_INVALID_INPUTS;
  }

  return ADJ_OK;
}

int adj_set_option(adj_adjointer* adjointer, int option, int choice)
{
  if (option < 0 || option >= ADJ_NO_OPTIONS)
  {
    strncpy(adj_error_msg, "Unknown option.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  adjointer->options[option] = choice;
  return ADJ_OK;
}

int adj_equation_count(adj_adjointer* adjointer, int* count)
{
  *count = adjointer->nequations;
  return ADJ_OK;
}

int adj_record_variable(adj_adjointer* adjointer, adj_variable var, adj_storage_data storage)
{
  adj_variable_data* data_ptr;
  int ierr;
  int compare_ierr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;

  ierr = adj_find_variable_data(&(adjointer->varhash), &var, &data_ptr);
  if (ierr != ADJ_OK && ierr != ADJ_ERR_HASH_FAILED) return ierr;

  if (ierr == ADJ_ERR_HASH_FAILED)
  {
    /* it's alright that this is the first time we've ever seen it */
    adj_variable_data* new_data;
    ierr = adj_add_new_hash_entry(adjointer, &var, &new_data);
    if (ierr != ADJ_OK) return ierr;
    new_data->equation = -1; /* it doesn't have an equation */
    data_ptr = new_data;
  }

  assert(data_ptr != NULL);

  if (storage.storage_memory_has_value && !data_ptr->storage.storage_memory_has_value) /* If we don't have a value recorded, any compare or overwrite flags can be ignored */
  {
  	return adj_record_variable_core_memory(adjointer, data_ptr, storage);
  }
  else if (storage.storage_disk_has_value && !data_ptr->storage.storage_disk_has_value) /* If we don't have a value recorded, any compare or overwrite flags can be ignored */
  {
    /* Generate the filename */
  	char filename[ADJ_NAME_LEN];
  	adj_variable_str(var, filename, ADJ_NAME_LEN);
  	strncat(filename, ".dat", 4);
  	strncpy(storage.storage_disk_filename, filename, ADJ_NAME_LEN);

    return adj_record_variable_core_disk(adjointer, data_ptr, storage);
  }
  else
  /* Sorry for the slight mess. The easiest way to understand this block is to build a 2x2 graph of
     compare and overwrite:

    |------------------|------------------|------------------|
    |                  |   Overwrite on   |  Overwrite off   |
    |------------------|------------------|------------------|
    |   Compare on     | Compare, then    | Just compare,    |
    |                  | overwrite; if the| and return the   |
    |                  | overwriting went | result of the    |
    |                  | OK, return the   | comparison       |
    |                  | result of the    |                  |
    |                  | comparison       |                  |
    |------------------|------------------|------------------|
    |   Compare off    | Just overwrite   | Return           |
    |                  |                  | ADJ_WARN_ALREA   |
    |                  |                  | DY_RECORDED      |
    |------------------|------------------|------------------| */
  {
    if (storage.compare)
    {
      if (!storage.storage_memory_has_value) /* Comparing variables only works when recording to memory */ 
      {
        strncpy(adj_error_msg, "Overwriting of values and comparing with previously computed values currently only works when recording to memory.", ADJ_ERROR_MSG_BUF);
        return ADJ_ERR_NOT_IMPLEMENTED;
      }

      ierr = adj_record_variable_compare(adjointer, data_ptr, var, storage);
      compare_ierr = ierr;
      if (storage.overwrite)
      {
        int record_ierr;

        ierr = adj_forget_variable_value(adjointer, data_ptr);
        if (ierr != ADJ_OK) return ierr;

        record_ierr = adj_record_variable_core_memory(adjointer, data_ptr, storage); /* Overwrite the result anyway */
        /* If no error happened from the recording, return the warning that the comparison failed;
           otherwise, return the (presumably more serious) error from the recording */
        if (record_ierr == ADJ_OK)
          return compare_ierr;
        else
          return record_ierr;
      }
      else /* We don't have the overwrite flag */
      {
        return ierr; /* Return the output of the comparison straight away */
      }
    }
    else /* We don't have the compare flag */
    {
      if (storage.overwrite)
      {
        ierr = adj_forget_variable_value(adjointer, data_ptr);
        if (ierr != ADJ_OK) return ierr;
        return adj_record_variable_core_memory(adjointer, data_ptr, storage);
      }
      else /* We don't have the overwrite flag */
      {
        char buf[ADJ_NAME_LEN];
        adj_variable_str(var, buf, ADJ_NAME_LEN);
        snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Variable %s already has a value.", buf);
        return ADJ_WARN_ALREADY_RECORDED;
      }
    }
  }

  return ADJ_OK; /* Should never get here, but keep the compiler quiet */
}

int adj_record_variable_core_memory(adj_adjointer* adjointer, adj_variable_data* data_ptr, adj_storage_data storage)
{
  if (adjointer->callbacks.vec_duplicate == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You have asked to record a value, but no ADJ_VEC_DUPLICATE_CB callback has been provided.");
    return ADJ_ERR_NEED_CALLBACK;
  }
  if (adjointer->callbacks.vec_axpy == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You have asked to record a value, but no ADJ_VEC_AXPY_CB callback has been provided.");
    return ADJ_ERR_NEED_CALLBACK;
  }

  assert(storage.storage_memory_has_value);

  switch (storage.storage_memory_type)
  {
    case ADJ_STORAGE_MEMORY_COPY:
      data_ptr->storage.storage_memory_type = ADJ_STORAGE_MEMORY_COPY;
      data_ptr->storage.storage_memory_has_value = storage.storage_memory_has_value;
      data_ptr->storage.storage_memory_is_checkpoint = storage.storage_memory_is_checkpoint;
      adjointer->callbacks.vec_duplicate(storage.value, &(data_ptr->storage.value));
      adjointer->callbacks.vec_axpy(&(data_ptr->storage.value), (adj_scalar)1.0, storage.value);
      break;
    case ADJ_STORAGE_MEMORY_INCREF:
      data_ptr->storage.storage_memory_type = ADJ_STORAGE_MEMORY_INCREF;
      data_ptr->storage.storage_memory_has_value = storage.storage_memory_has_value;
      data_ptr->storage.storage_memory_is_checkpoint = storage.storage_memory_is_checkpoint;
      break;
    default:
      strncpy(adj_error_msg, "Memory storage types other than ADJ_STORAGE_MEMORY_COPY and ADJ_STORAGE_MEMORY_INCREF  are not implemented yet.", ADJ_ERROR_MSG_BUF);
      return ADJ_ERR_NOT_IMPLEMENTED;
  }

  return ADJ_OK;
}

int adj_record_variable_core_disk(adj_adjointer* adjointer, adj_variable_data* data_ptr, adj_storage_data storage)
{
  if (adjointer->callbacks.vec_to_file == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You have asked to record a value to disk, but no ADJ_VEC_TO_FILE_CB callback has been provided.");
    return ADJ_ERR_NEED_CALLBACK;
  }

  assert(storage.storage_disk_has_value);

  data_ptr->storage.storage_disk_has_value = storage.storage_disk_has_value;
  data_ptr->storage.storage_disk_is_checkpoint = storage.storage_disk_is_checkpoint;
  strncpy(data_ptr->storage.storage_disk_filename, storage.storage_disk_filename, ADJ_NAME_LEN);
  adjointer->callbacks.vec_to_file(storage.value, data_ptr->storage.storage_disk_filename);

  return ADJ_OK;
}

int adj_record_variable_compare(adj_adjointer* adjointer, adj_variable_data* data_ptr, adj_variable var, adj_storage_data storage)
{
  if (data_ptr->storage.storage_memory_has_value && storage.compare)
  {
    adj_vector tmp;
    adj_scalar norm;

    if (adjointer->callbacks.vec_get_norm == NULL)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You have asked to compare a value against one already recorded, but no ADJ_VEC_GET_NORM_CB callback has been provided.");
      return ADJ_ERR_NEED_CALLBACK;
    }
    if (adjointer->callbacks.vec_duplicate == NULL)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You have asked to compare a value against one already recorded, but no ADJ_VEC_DUPLICATE_CB callback has been provided.");
      return ADJ_ERR_NEED_CALLBACK;
    }
    if (adjointer->callbacks.vec_axpy == NULL)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You have asked to compare a value against one already recorded, but no ADJ_VEC_AXPY_CB callback has been provided.");
      return ADJ_ERR_NEED_CALLBACK;
    }
    if (adjointer->callbacks.vec_destroy == NULL)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You have asked to compare a value against one already recorded, but no ADJ_VEC_DESTROY_CB callback has been provided.");
      return ADJ_ERR_NEED_CALLBACK;
    }
    if (storage.storage_memory_type != ADJ_STORAGE_MEMORY_COPY && storage.storage_memory_type != ADJ_STORAGE_MEMORY_INCREF)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Sorry, comparison of values hasn't been generalised to storage types other than ADJ_STORAGE_MEMORY_COPY and ADJ_STORAGE_MEMORY_INCREF yet.");
      return ADJ_ERR_NOT_IMPLEMENTED;
      /* For future developers: the reason is that storage.value (used a few lines below) might not exist */
    }

    adjointer->callbacks.vec_duplicate(data_ptr->storage.value, &tmp);
    adjointer->callbacks.vec_axpy(&tmp, (adj_scalar)1.0, data_ptr->storage.value);
    adjointer->callbacks.vec_axpy(&tmp, (adj_scalar)-1.0, storage.value);
    adjointer->callbacks.vec_get_norm(tmp, &norm);
    adjointer->callbacks.vec_destroy(&tmp);

    if (norm > storage.comparison_tolerance) /* Greater than, so that we can use a comparison tolerance of 0.0 */
    {
      char buf[ADJ_NAME_LEN];
      adj_variable_str(var, buf, ADJ_NAME_LEN);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Comparing %s against previously recorded value: norm of the difference is %e (> tolerance of %e)", buf, norm, storage.comparison_tolerance);
      return ADJ_WARN_COMPARISON_FAILED;
    }
    else
    {
      return ADJ_OK;
    }
  }
  else
  {
    return ADJ_OK;
  }
}

int adj_register_operator_callback(adj_adjointer* adjointer, int type, char* name, void (*fn)(void))
{
  adj_op_callback_list* cb_list_ptr;
  adj_op_callback* cb_ptr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;

  switch(type)
  {
    case ADJ_NBLOCK_COLOURING_CB:
      cb_list_ptr = &(adjointer->nonlinear_colouring_list);
      break;
    case ADJ_NBLOCK_ACTION_CB:
      cb_list_ptr = &(adjointer->nonlinear_action_list);
      break;
    case ADJ_NBLOCK_DERIVATIVE_ACTION_CB:
      cb_list_ptr = &(adjointer->nonlinear_derivative_action_list);
      break;
    case ADJ_NBLOCK_DERIVATIVE_ASSEMBLY_CB:
      cb_list_ptr = &(adjointer->nonlinear_derivative_assembly_list);
      break;
    case ADJ_BLOCK_ACTION_CB:
      cb_list_ptr = &(adjointer->block_action_list);
      break;
    case ADJ_BLOCK_ASSEMBLY_CB:
      cb_list_ptr = &(adjointer->block_assembly_list);
      break;
    default:
      strncpy(adj_error_msg, "Unknown callback type.", ADJ_ERROR_MSG_BUF);
      return ADJ_ERR_INVALID_INPUTS;
  }
  /* First, we look for an existing callback data structure that might already exist, to replace the function */
  cb_ptr = cb_list_ptr->firstnode;
  while (cb_ptr != NULL)
  {
    if (strncmp(cb_ptr->name, name, ADJ_NAME_LEN) == 0)
    {
      cb_ptr->callback = fn;
      return ADJ_OK;
    }
    cb_ptr = cb_ptr->next;
  }

  /* If we got here, that means that we didn't find it. Tack it on to the end of the list. */
  cb_ptr = (adj_op_callback*) malloc(sizeof(adj_op_callback));
  ADJ_CHKMALLOC(cb_ptr);
  strncpy(cb_ptr->name, name, ADJ_NAME_LEN);
  cb_ptr->callback = fn;
  cb_ptr->next = NULL;

  /* Special case for the first callback */
  if (cb_list_ptr->firstnode == NULL)
  {
    cb_list_ptr->firstnode = cb_ptr;
    cb_list_ptr->lastnode = cb_ptr;
  }
  else
  {
    cb_list_ptr->lastnode->next = cb_ptr;
    cb_list_ptr->lastnode = cb_ptr;
  }

  return ADJ_OK;
}

int adj_register_data_callback(adj_adjointer* adjointer, int type, void (*fn)(void))
{
  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;

  switch (type)
  {
    case ADJ_VEC_DUPLICATE_CB:
      adjointer->callbacks.vec_duplicate = (void(*)(adj_vector x, adj_vector *newx)) fn;
      break;
    case ADJ_VEC_AXPY_CB:
      adjointer->callbacks.vec_axpy = (void(*)(adj_vector *y, adj_scalar alpha, adj_vector x)) fn;
      break;
    case ADJ_VEC_DESTROY_CB:
      adjointer->callbacks.vec_destroy = (void(*)(adj_vector*)) fn;
      break;
    case ADJ_VEC_DIVIDE_CB:
      adjointer->callbacks.vec_divide = (void(*)(adj_vector *numerator, adj_vector denominator)) fn;
      break;
    case ADJ_VEC_SET_VALUES_CB:
      adjointer->callbacks.vec_set_values = (void(*)(adj_vector *vec, adj_scalar scalars[])) fn;
      break;
    case ADJ_VEC_GET_SIZE_CB:
      adjointer->callbacks.vec_get_size = (void(*)(adj_vector vec, int *sz)) fn;
      break;
    case ADJ_VEC_GET_NORM_CB:
      adjointer->callbacks.vec_get_norm = (void(*)(adj_vector vec, adj_scalar *norm)) fn;
      break;
    case ADJ_VEC_DOT_PRODUCT_CB:
      adjointer->callbacks.vec_dot_product = (void(*)(adj_vector x, adj_vector y, adj_scalar* val)) fn;
      break;
    case ADJ_VEC_SET_RANDOM_CB:
      adjointer->callbacks.vec_set_random = (void(*)(adj_vector* x)) fn;
      break;
    case ADJ_VEC_TO_FILE_CB:
      adjointer->callbacks.vec_to_file = (void(*)(adj_vector x, char* filename)) fn;
      break;
    case ADJ_VEC_FROM_FILE_CB:
      adjointer->callbacks.vec_from_file = (void(*)(adj_vector* x, char* filename)) fn;
      break;

    case ADJ_MAT_DUPLICATE_CB:
      adjointer->callbacks.mat_duplicate = (void(*)(adj_matrix matin, adj_matrix *matout)) fn;
      break;
    case ADJ_MAT_AXPY_CB:
      adjointer->callbacks.mat_axpy = (void(*)(adj_matrix *Y, adj_scalar alpha, adj_matrix X)) fn;
      break;
    case ADJ_MAT_DESTROY_CB:
      adjointer->callbacks.mat_destroy = (void(*)(adj_matrix *mat)) fn;
      break;
    case ADJ_SOLVE_CB:
      adjointer->callbacks.solve = (void(*)(adj_variable var, adj_matrix mat, adj_vector rhs, adj_vector *soln)) fn;
      break;

   default:
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Unknown data callback type %d.", type);
      return ADJ_ERR_INVALID_INPUTS;
  }

  return ADJ_OK;
}

int adj_register_functional_callback(adj_adjointer* adjointer, char* name, void (*fn)(adj_adjointer* adjointer, int timestep, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output))
{
  adj_func_callback_list* cb_list_ptr;
  adj_func_callback* cb_ptr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;

  cb_list_ptr = &(adjointer->functional_list);

  /* First, we look for an existing callback data structure that might already exist, to replace the function */
  cb_ptr = cb_list_ptr->firstnode;
  while (cb_ptr != NULL)
  {
    if (strncmp(cb_ptr->name, name, ADJ_NAME_LEN) == 0)
    {
      cb_ptr->callback = (void (*)(void* adjointer, int timestep, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output)) fn;
      return ADJ_OK;
    }
    cb_ptr = cb_ptr->next;
  }

  /* If we got here, that means that we didn't find it. Tack it on to the end of the list. */
  cb_ptr = (adj_func_callback*) malloc(sizeof(adj_func_callback));
  ADJ_CHKMALLOC(cb_ptr);
  strncpy(cb_ptr->name, name, ADJ_NAME_LEN);
  cb_ptr->callback = (void (*)(void* adjointer, int timestep, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output)) fn;
  cb_ptr->next = NULL;

  /* Special case for the first callback */
  if (cb_list_ptr->firstnode == NULL)
  {
    cb_list_ptr->firstnode = cb_ptr;
    cb_list_ptr->lastnode = cb_ptr;
  }
  else
  {
    cb_list_ptr->lastnode->next = cb_ptr;
    cb_list_ptr->lastnode = cb_ptr;
  }

  return ADJ_OK;
}

int adj_register_functional_derivative_callback(adj_adjointer* adjointer, char* name, void (*fn)(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output))
{
  adj_func_deriv_callback_list* cb_list_ptr;
  adj_func_deriv_callback* cb_ptr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;

  cb_list_ptr = &(adjointer->functional_derivative_list);

  /* First, we look for an existing callback data structure that might already exist, to replace the function */
  cb_ptr = cb_list_ptr->firstnode;
  while (cb_ptr != NULL)
  {
    if (strncmp(cb_ptr->name, name, ADJ_NAME_LEN) == 0)
    {
      cb_ptr->callback = (void (*)(void* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output)) fn;
      return ADJ_OK;
    }
    cb_ptr = cb_ptr->next;
  }

  /* If we got here, that means that we didn't find it. Tack it on to the end of the list. */
  cb_ptr = (adj_func_deriv_callback*) malloc(sizeof(adj_func_deriv_callback));
  ADJ_CHKMALLOC(cb_ptr);
  strncpy(cb_ptr->name, name, ADJ_NAME_LEN);
  cb_ptr->callback = (void (*)(void* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output)) fn;
  cb_ptr->next = NULL;

  /* Special case for the first callback */
  if (cb_list_ptr->firstnode == NULL)
  {
    cb_list_ptr->firstnode = cb_ptr;
    cb_list_ptr->lastnode = cb_ptr;
  }
  else
  {
    cb_list_ptr->lastnode->next = cb_ptr;
    cb_list_ptr->lastnode = cb_ptr;
  }

  return ADJ_OK;
}

int adj_register_forward_source_callback(adj_adjointer* adjointer, void (*fn)(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, void* context, adj_vector* output, int* has_output))
{
  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;
  adjointer->forward_source_callback = fn;
  return ADJ_OK;
}

int adj_forget_adjoint_equation(adj_adjointer* adjointer, int equation)
{
  adj_variable_data* data;
  int should_we_delete;
  int i;
  int ierr;
  int min_timestep;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;

  if (equation >= adjointer->nequations)
  {
    strncpy(adj_error_msg, "No such equation.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  data = adjointer->vardata.firstnode;
  while (data != NULL)
  {
    if (data->storage.storage_memory_has_value || data->storage.storage_disk_has_value)
    {
      should_we_delete = 1;
      /* Check the adjoint equations we could explicitly compute */
      for (i = 0; i < data->nadjoint_equations; i++)
      {
        if (equation > data->adjoint_equations[i])
        {
          should_we_delete = 0;
          break;
        }
      }

      /* Also check that it isn't necessary for any timesteps we still have to compute
         functional right-hand-sides for */
      if (data->ndepending_timesteps > 0 && should_we_delete == 1)
      {
        int start_equation;
        min_timestep = adj_minval(data->depending_timesteps, data->ndepending_timesteps);
        ierr = adj_timestep_start_equation(adjointer, min_timestep, &start_equation);
        assert(ierr == ADJ_OK);

        if (equation > start_equation)
        {
          should_we_delete = 0;
        }
      }

      if (should_we_delete)
      {
        ierr = adj_forget_variable_value(adjointer, data);
        if (ierr != ADJ_OK) return ierr;
      }
    }

    data = data->next;
  }

  return ADJ_OK;
}

/*
 * Forgets the forward variables that are not needed for any forward equations larger than "equation"
 */
int adj_forget_forward_equation(adj_adjointer* adjointer, int equation)
{
  int last_equation, ierr;
  ierr = adj_equation_count(adjointer, &last_equation);
  if (ierr != ADJ_OK) return ierr;

	return adj_forget_forward_equation_until(adjointer, equation, last_equation);
}

/*
 * Like adj_forget_forward_equation, but assumes that the annotation stops at last_equation.
 */
int adj_forget_forward_equation_until(adj_adjointer* adjointer, int equation, int last_equation)
{
  adj_variable_data* data;
  int should_we_delete;
  int i;
  int ierr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;

  if (equation >= adjointer->nequations)
  {
    strncpy(adj_error_msg, "No such equation.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  data = adjointer->vardata.firstnode;
  while (data != NULL)
  {
  	/* Only forget forward variables */
  	if (data->equation<0) /* Adjoint variables have no forward equation */
  	{
  		data = data->next;
  		continue;
  	}

    if (data->storage.storage_memory_has_value || data->storage.storage_disk_has_value)
    {
      should_we_delete = 1;
      /* Check the forward equations we could explicitly compute */
      for (i = 0; i < data->ntargeting_equations; i++)
      {
        if (equation < data->targeting_equations[i] && data->targeting_equations[i] <= last_equation)
        {
          should_we_delete = 0;
          break;
        }
      }

      for (i = 0; i < data->ndepending_equations; i++)
      {
        if (equation < data->depending_equations[i] && data->depending_equations[i] <= last_equation)
        {
          should_we_delete = 0;
          break;
        }
      }

      /* Also check that it isn't necessary for any timesteps we still have to compute
         right-hand-sides for */
      for (i = 0; i < data->nrhs_equations; i++)
      {
        if (equation < data->rhs_equations[i] && data->rhs_equations[i] <= last_equation)
        {
          should_we_delete = 0;
          break;
        }
      }

      if (should_we_delete)
      {
      	/* Forget variables that are not checkpoints */
      	if (data->storage.storage_disk_has_value && !data->storage.storage_disk_is_checkpoint)
      	{
      		ierr = adj_forget_variable_value_from_disk(adjointer, data);
      		if (ierr != ADJ_OK) return ierr;
      	}
      	if (data->storage.storage_memory_has_value && !data->storage.storage_memory_is_checkpoint)
      	{
      		ierr = adj_forget_variable_value_from_memory(adjointer, data);
      		if (ierr != ADJ_OK) return ierr;
      	}
      }
    }

    data = data->next;
  }

  return ADJ_OK;
}
int adj_find_operator_callback(adj_adjointer* adjointer, int type, char* name, void (**fn)(void))
{
  adj_op_callback_list* cb_list_ptr;
  adj_op_callback* cb_ptr;

  char adj_callback_types[6][ADJ_ERROR_MSG_BUF] = {"ADJ_NBLOCK_COLOURING_CB", "ADJ_NBLOCK_ACTION_CB", "ADJ_NBLOCK_DERIVATIVE_ACTION_CB",
                                                   "ADJ_NBLOCK_DERIVATIVE_ASSEMBLY_CB", "ADJ_BLOCK_ACTION_CB", "ADJ_BLOCK_ASSEMBLY_CB"};

  switch(type)
  {
    case ADJ_NBLOCK_COLOURING_CB:
      cb_list_ptr = &(adjointer->nonlinear_colouring_list);
      break;
    case ADJ_NBLOCK_ACTION_CB:
      cb_list_ptr = &(adjointer->nonlinear_action_list);
      break;
    case ADJ_NBLOCK_DERIVATIVE_ACTION_CB:
      cb_list_ptr = &(adjointer->nonlinear_derivative_action_list);
      break;
    case ADJ_NBLOCK_DERIVATIVE_ASSEMBLY_CB:
      cb_list_ptr = &(adjointer->nonlinear_derivative_assembly_list);
      break;
    case ADJ_BLOCK_ACTION_CB:
      cb_list_ptr = &(adjointer->block_action_list);
      break;
    case ADJ_BLOCK_ASSEMBLY_CB:
      cb_list_ptr = &(adjointer->block_assembly_list);
      break;
    default:
      strncpy(adj_error_msg, "Unknown callback type.", ADJ_ERROR_MSG_BUF);
      return ADJ_ERR_INVALID_INPUTS;
  }

  cb_ptr = cb_list_ptr->firstnode;
  while (cb_ptr != NULL)
  {
    if (strncmp(cb_ptr->name, name, ADJ_NAME_LEN) == 0)
    {
      *fn = cb_ptr->callback;
      return ADJ_OK;
    }
    cb_ptr = cb_ptr->next;
  }

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Could not find callback %s for operator %s.", adj_callback_types[type-1], name);
  return ADJ_ERR_NEED_CALLBACK;
}

int adj_find_functional_callback(adj_adjointer* adjointer, char* name, void (**fn)(adj_adjointer* adjointer, int timestep, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output))
{
  adj_func_callback_list* cb_list_ptr;
  adj_func_callback* cb_ptr;

  cb_list_ptr = &(adjointer->functional_list);

  cb_ptr = cb_list_ptr->firstnode;
  while (cb_ptr != NULL)
  {
    if (strncmp(cb_ptr->name, name, ADJ_NAME_LEN) == 0)
    {
      *fn = (void (*)(adj_adjointer* adjointer, int timestep, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output)) cb_ptr->callback;
      return ADJ_OK;
    }
    cb_ptr = cb_ptr->next;
  }

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Could not find functional callback %s.", name);
  return ADJ_ERR_NEED_CALLBACK;
}

int adj_find_functional_derivative_callback(adj_adjointer* adjointer, char* name, void (**fn)(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output))
{
  adj_func_deriv_callback_list* cb_list_ptr;
  adj_func_deriv_callback* cb_ptr;

  cb_list_ptr = &(adjointer->functional_derivative_list);

  cb_ptr = cb_list_ptr->firstnode;
  while (cb_ptr != NULL)
  {
    if (strncmp(cb_ptr->name, name, ADJ_NAME_LEN) == 0)
    {
      *fn = (void (*)(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output)) cb_ptr->callback;
      return ADJ_OK;
    }
    cb_ptr = cb_ptr->next;
  }

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Could not find functional derivative callback %s.", name);
  return ADJ_ERR_NEED_CALLBACK;
}

int adj_get_variable_value(adj_adjointer* adjointer, adj_variable var, adj_vector* value)
{
  int ierr;
  adj_variable_data* data_ptr;

  ierr = adj_find_variable_data(&(adjointer->varhash), &var, &data_ptr);
  if (ierr != ADJ_OK) return ierr;

  if (!data_ptr->storage.storage_memory_has_value && !data_ptr->storage.storage_disk_has_value)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);

    ierr = ADJ_ERR_NEED_VALUE;
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but don't have one recorded.", buf);
    return ierr;
  }

  /* Memory storage */
  if (data_ptr->storage.storage_memory_has_value)
  {
   *value = data_ptr->storage.value;
   return ADJ_OK;
  }
  /* Disk storage */
  else if(data_ptr->storage.storage_disk_has_value)
  {
    if (adjointer->callbacks.vec_from_file == NULL)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You have asked to get a value from disk, but no ADJ_VEC_FROM_FILE_CB callback has been provided.");
      return ADJ_ERR_NEED_CALLBACK;
    }
    adjointer->callbacks.vec_from_file(value, data_ptr->storage.storage_disk_filename);
  }
  return ADJ_OK;
}

int adj_has_variable_value(adj_adjointer* adjointer, adj_variable var)
{
  if((adj_has_variable_value_memory(adjointer, var)==ADJ_OK) || (adj_has_variable_value_disk(adjointer, var)==ADJ_OK))
  	return ADJ_OK;
  else
  	return ADJ_ERR_NEED_VALUE;
}

int adj_has_variable_value_memory(adj_adjointer* adjointer, adj_variable var)
{
  int ierr;
  adj_variable_data* data_ptr;

  ierr = adj_find_variable_data(&(adjointer->varhash), &var, &data_ptr);
  if (ierr != ADJ_OK)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but have never seen it.", buf);
    return ierr;
  }

  if (!data_ptr->storage.storage_memory_has_value)
    return ADJ_ERR_NEED_VALUE;

  return ADJ_OK;
}

int adj_has_variable_value_disk(adj_adjointer* adjointer, adj_variable var)
{
  int ierr;
  adj_variable_data* data_ptr;

  ierr = adj_find_variable_data(&(adjointer->varhash), &var, &data_ptr);
  if (ierr != ADJ_OK)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but have never seen it.", buf);
    return ierr;
  }

  if (!data_ptr->storage.storage_disk_has_value)
    return ADJ_ERR_NEED_VALUE;

  return ADJ_OK;
}

int adj_is_checkpoint_variable_disk(adj_adjointer* adjointer, adj_variable var)
{
  int ierr;
  adj_variable_data* data_ptr;

  ierr = adj_find_variable_data(&(adjointer->varhash), &var, &data_ptr);
  if (ierr != ADJ_OK)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but have never seen it.", buf);
    return ierr;
  }

  return data_ptr->storage.storage_disk_is_checkpoint;
}

int adj_is_checkpoint_variable_memory(adj_adjointer* adjointer, adj_variable var)
{
  int ierr;
  adj_variable_data* data_ptr;

  ierr = adj_find_variable_data(&(adjointer->varhash), &var, &data_ptr);
  if (ierr != ADJ_OK)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but have never seen it.", buf);
    return ierr;
  }

  return data_ptr->storage.storage_memory_is_checkpoint;
}

int adj_forget_variable_value(adj_adjointer* adjointer, adj_variable_data* data)
{
  int ierr;

  if (data->storage.storage_disk_has_value)
  {
    ierr = adj_forget_variable_value_from_disk(adjointer, data);
    if (ierr != ADJ_OK) return ierr;
  }
  if (data->storage.storage_memory_has_value)
  {
    ierr = adj_forget_variable_value_from_memory(adjointer, data);
    if (ierr != ADJ_OK) return ierr;
  }
  return ADJ_OK;
}

int adj_forget_variable_value_from_disk(adj_adjointer* adjointer, adj_variable_data* data)
{
  FILE *istream;
  int ierr;

  assert(data->storage.storage_disk_has_value);
  if (!(istream = fopen(data->storage.storage_disk_filename, "r+"))) { 
    char buf[ADJ_NAME_LEN];
    adj_variable_str(adjointer->equations[data->equation].variable, buf, ADJ_NAME_LEN);

    ierr = ADJ_ERR_INVALID_INPUTS;
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Can not access variable %s in file '%s'.", buf, data->storage.storage_disk_filename);
    return ierr;
  }
  fclose(istream); 

  data->storage.storage_disk_has_value = ADJ_FALSE;
  remove(data->storage.storage_disk_filename);
  return ADJ_OK;
}

int adj_forget_variable_value_from_memory(adj_adjointer* adjointer, adj_variable_data* data)
{
  if (adjointer->callbacks.vec_destroy == NULL)
  {
    strncpy(adj_error_msg, "Need vec_destroy data callback.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_NEED_CALLBACK;
  }

  assert(data->storage.storage_memory_has_value);

  data->storage.storage_memory_has_value = ADJ_FALSE;
  adjointer->callbacks.vec_destroy(&(data->storage.value));
  return ADJ_OK;
}

int adj_destroy_variable_data(adj_adjointer* adjointer, adj_variable_data* data)
{

  if (data->ntargeting_equations > 0)
  {
    free(data->targeting_equations);
    data->ntargeting_equations = 0;
  }

  if (data->ndepending_equations > 0)
  {
    free(data->depending_equations);
    data->ndepending_equations = 0;
  }

  if (data->nrhs_equations > 0)
  {
    free(data->rhs_equations);
    data->nrhs_equations = 0;
  }

  if (data->nadjoint_equations > 0)
  {
    free(data->adjoint_equations);
    data->nadjoint_equations = 0;
  }

  if (data->storage.storage_memory_has_value || data->storage.storage_disk_has_value)
    return adj_forget_variable_value(adjointer, data);

  return ADJ_OK;
}

int adj_storage_memory_copy(adj_vector value, adj_storage_data* data)
{
  data->storage_memory_has_value = ADJ_TRUE;
  data->storage_memory_type = ADJ_STORAGE_MEMORY_COPY;
  data->value = value;

  data->storage_disk_has_value = ADJ_FALSE;

  data->compare = ADJ_FALSE;
  data->comparison_tolerance = (adj_scalar)0.0;
  data->overwrite = ADJ_FALSE;
  data->storage_memory_is_checkpoint = ADJ_FALSE;
  data->storage_disk_is_checkpoint = ADJ_FALSE;
  return ADJ_OK;
}

int adj_storage_memory_incref(adj_vector value, adj_storage_data* data)
{
  data->storage_memory_has_value = ADJ_TRUE;
  data->storage_memory_type = ADJ_STORAGE_MEMORY_INCREF;
  data->value = value;

  data->storage_disk_has_value = ADJ_FALSE;

  data->compare = ADJ_FALSE;
  data->comparison_tolerance = (adj_scalar)0.0;
  data->overwrite = ADJ_FALSE;
  data->storage_memory_is_checkpoint = ADJ_FALSE;
  data->storage_disk_is_checkpoint = ADJ_FALSE;
  return ADJ_OK;
}

int adj_storage_disk(adj_vector value, adj_storage_data* data)
{
  data->value = value;
  data->storage_memory_has_value = ADJ_FALSE;
  data->storage_memory_type = ADJ_UNSET;

  data->storage_disk_has_value = ADJ_TRUE;

  data->compare = ADJ_FALSE;
  data->comparison_tolerance = (adj_scalar)0.0;
  data->overwrite = ADJ_FALSE;
  data->storage_memory_is_checkpoint = ADJ_FALSE;
  data->storage_disk_is_checkpoint = ADJ_FALSE;
  return ADJ_OK;
}

int adj_storage_set_compare(adj_storage_data* data, int compare, adj_scalar comparison_tolerance)
{
  if (compare != ADJ_TRUE && compare != ADJ_FALSE)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "compare must either be ADJ_TRUE or ADJ_FALSE.");
    return ADJ_ERR_INVALID_INPUTS;
  }

  if ((float) comparison_tolerance < 0.0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "comparison_tolerance must be >= 0.0.");
    return ADJ_ERR_INVALID_INPUTS;
  }

  data->compare = compare;
  data->comparison_tolerance = comparison_tolerance;
  return ADJ_OK;
}

int adj_storage_set_overwrite(adj_storage_data* data, int overwrite)
{
  if (overwrite != ADJ_TRUE && overwrite != ADJ_FALSE)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "overwrite must either be ADJ_TRUE or ADJ_FALSE.");
    return ADJ_ERR_INVALID_INPUTS;
  }

  data->overwrite = overwrite;
  return ADJ_OK;
}

int adj_storage_set_checkpoint(adj_storage_data* data, int checkpoint)
{
	if (checkpoint != ADJ_TRUE && checkpoint != ADJ_FALSE)
	{
	  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "checkpoint must either be ADJ_TRUE or ADJ_FALSE.");
	  return ADJ_ERR_INVALID_INPUTS;
	}

	data->storage_memory_is_checkpoint=checkpoint;
	data->storage_disk_is_checkpoint=checkpoint;
	return ADJ_OK;
}

int adj_add_new_hash_entry(adj_adjointer* adjointer, adj_variable* var, adj_variable_data** data)
{
  int ierr;
  adj_variable auxvar;
  adj_variable_data* aux_data;

  /* First, check that we don't have a corresponding variable with the opposite sense
     of auxiliary-ness already registered -- if so, the user has probably made a mistake */
  auxvar = *var;
  auxvar.auxiliary = !auxvar.auxiliary;
  ierr = adj_find_variable_data(&(adjointer->varhash), &auxvar, &aux_data);
  if (ierr != ADJ_ERR_HASH_FAILED)
  {
    char buf[ADJ_NAME_LEN];
    char auxbuf[ADJ_NAME_LEN];
    adj_variable_str(*var, buf, ADJ_NAME_LEN);
    adj_variable_str(auxvar, auxbuf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Tried to add a hash entry for %s, but already have a hash entry for %s.", buf, auxbuf);
    return ADJ_ERR_INVALID_INPUTS;
  }

  *data = (adj_variable_data*) malloc(sizeof(adj_variable_data));
  ADJ_CHKMALLOC(*data);
  memset(*data, 0, sizeof(adj_variable_data));
  (*data)->equation = -1;
  (*data)->next = NULL;
  (*data)->storage.storage_memory_has_value = 0;
  (*data)->storage.storage_memory_type = ADJ_UNSET;
  (*data)->storage.storage_disk_has_value = 0;
  (*data)->ntargeting_equations = 0;
  (*data)->targeting_equations = NULL;
  (*data)->ndepending_equations = 0;
  (*data)->depending_equations = NULL;
  (*data)->nrhs_equations = 0;
  (*data)->rhs_equations = NULL;
  (*data)->ndepending_timesteps = 0;
  (*data)->depending_timesteps = NULL;
  (*data)->nadjoint_equations = 0;
  (*data)->adjoint_equations = NULL;

  /* add to the hash table */
  ierr = adj_add_variable_data(&(adjointer->varhash), var, *data);
  if (ierr != ADJ_OK) return ierr;

  /* and add to the data list */
  if (adjointer->vardata.firstnode == NULL)
  {
    adjointer->vardata.firstnode = *data;
    adjointer->vardata.lastnode = *data;
  }
  else
  {
    adjointer->vardata.lastnode->next = *data;
    adjointer->vardata.lastnode = *data;
  }

  return ADJ_OK;
}

int adj_timestep_count(adj_adjointer* adjointer, int* count)
{
  *count = adjointer->ntimesteps;
  return ADJ_OK;
}

int adj_iteration_count(adj_adjointer* adjointer, adj_variable variable, int* count)
{
  int ierr;
  adj_variable_data* data;
  *count = -1;
  do
  {
    *count=*count+1;
    variable.iteration = *count;
    ierr = adj_find_variable_data(&(adjointer->varhash), &variable, &data);
  }
  while (ierr == ADJ_OK);

  if (*count==0)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(variable, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Error in adj_iteration_count: No iteration found for supplied variable %s.", buf);
    return ADJ_ERR_INVALID_INPUTS;
  }
  return ADJ_OK;
}

int adj_timestep_start_equation(adj_adjointer* adjointer, int timestep, int* start)
{
  if (timestep < 0 || timestep >= adjointer->ntimesteps)
  {
    strncpy(adj_error_msg, "Invalid timestep supplied to adj_timestep_start.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  *start = adjointer->timestep_data[timestep].start_equation;
  return ADJ_OK;
}

int adj_timestep_end_equation(adj_adjointer* adjointer, int timestep, int* end)
{
  if (timestep < 0 || timestep >= adjointer->ntimesteps)
  {
    strncpy(adj_error_msg, "Invalid timestep supplied to adj_timestep_end.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  if (timestep < adjointer->ntimesteps-1)
  {
    *end = adjointer->timestep_data[timestep+1].start_equation - 1;
  }
  else
  {
    *end = adjointer->nequations - 1;
  }
  return ADJ_OK;
}

int adj_timestep_set_times(adj_adjointer* adjointer, int timestep, adj_scalar start, adj_scalar end)
{
  int ierr;

  if (timestep < 0)
  {
    strncpy(adj_error_msg, "Invalid timestep supplied to adj_timestep_set_times.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  if (end <= start)
  {
    strncpy(adj_error_msg, "End time cannot be less than start time.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  if (adjointer->ntimesteps <= timestep)
  {
    ierr = adj_extend_timestep_data(adjointer, timestep + 1);
    if (ierr != ADJ_OK) return ierr;
  }

  adjointer->timestep_data[timestep].start_time = start;
  adjointer->timestep_data[timestep].end_time = end;

  return ADJ_OK;
}

int adj_timestep_get_times(adj_adjointer* adjointer, int timestep, adj_scalar* start, adj_scalar* end)
{
  if (timestep < 0 || timestep >= adjointer->ntimesteps)
  {
    strncpy(adj_error_msg, "Invalid timestep supplied to adj_timestep_get_times.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  *start = adjointer->timestep_data[timestep].start_time;
  *end   = adjointer->timestep_data[timestep].end_time;

  /* A special exception for the last timestep, as it may not be a real
     timestep; it might only be introduced internally to be a container
     for the very last equation. In that case, we want to give it a sensible
     time, so that the user is not confused. */
  if (*start == ADJ_UNSET && *end == ADJ_UNSET 
      && timestep == adjointer->ntimesteps - 1
      && timestep - 1 >= 0)
  {
    *start = adjointer->timestep_data[timestep-1].end_time;
    *end   = adjointer->timestep_data[timestep-1].end_time;
  }

  if (*start == ADJ_UNSET && *end == ADJ_UNSET)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You have asked for the times for timestep %d, but they have not been set.", timestep);
    return ADJ_ERR_INVALID_INPUTS;
  }
  else
    return ADJ_OK;
}

int adj_timestep_set_functional_dependencies(adj_adjointer* adjointer, int timestep, char* functional, int ndepends, adj_variable* dependencies)
{
  int i;
  int ierr;

  adj_functional_data* functional_data_ptr = NULL;
  if (timestep < 0)
  {
    strncpy(adj_error_msg, "Invalid timestep supplied to adj_timestep_set_times.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  if (ndepends < 0)
  {
    strncpy(adj_error_msg, "Must supply a nonnegative number of dependencies.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  for (i = 0; i < ndepends; i++)
  {
    if (dependencies[i].type != ADJ_FORWARD)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Functional dependencies must be forward variables.");
      return ADJ_ERR_INVALID_INPUTS;
    }
  }

  if (adjointer->ntimesteps <= timestep)
  {
    ierr = adj_extend_timestep_data(adjointer, timestep + 1);
    if (ierr != ADJ_OK) return ierr;
  }

  /* Make sure that the dependencies for this timestep have not been set before */
  functional_data_ptr = adjointer->timestep_data[timestep].functional_data_start;
  while (functional_data_ptr != NULL)
  {
    if (strncmp(functional_data_ptr->name, functional, ADJ_NAME_LEN) == 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The dependencies for functional %s at timestep %d have been already set.", functional, timestep);
      return ADJ_ERR_INVALID_INPUTS;
    }
    functional_data_ptr = functional_data_ptr->next;
  }

  /* append the functional information to the adjointer */  
  functional_data_ptr = (adj_functional_data*) malloc(sizeof(adj_functional_data));
  ADJ_CHKMALLOC(functional_data_ptr);

  if (adjointer->timestep_data[timestep].functional_data_start == NULL)
    adjointer->timestep_data[timestep].functional_data_start = functional_data_ptr;
  if (adjointer->timestep_data[timestep].functional_data_end != NULL) 
    adjointer->timestep_data[timestep].functional_data_end->next = functional_data_ptr;
  adjointer->timestep_data[timestep].functional_data_end = functional_data_ptr;

  functional_data_ptr->next = NULL;
  strncpy(functional_data_ptr->name, functional, ADJ_NAME_LEN);
  functional_data_ptr->ndepends = ndepends;
  functional_data_ptr->dependencies = (adj_variable*) malloc(ndepends * sizeof(adj_variable));
  ADJ_CHKMALLOC(functional_data_ptr->dependencies);
  memcpy(functional_data_ptr->dependencies, dependencies, ndepends * sizeof(adj_variable));

  for (i = 0; i < ndepends; i++)
  {
    int ierr;
    adj_variable_data* data_ptr;

    ierr = adj_find_variable_data(&(adjointer->varhash), &(dependencies[i]), &data_ptr);
    if (ierr == ADJ_ERR_HASH_FAILED)
    {
      ierr = adj_add_new_hash_entry(adjointer, &(dependencies[i]), &data_ptr);
      if (ierr != ADJ_OK) return ierr;
    }

    /* Record that this variable is necessary for the functional evaluation at this point in time */
    ierr = adj_append_unique(&(data_ptr->depending_timesteps), &(data_ptr->ndepending_timesteps), timestep);
    if (ierr != ADJ_OK) return ierr;
  }
  /* We are done */
  return ADJ_OK;
}

int adj_append_unique(int** array, int* array_sz, int value)
{
  int i;

  for (i = 0; i < *array_sz; i++)
    if ((*array)[i] == value)
      return ADJ_OK;

  /* So if we got here, we really do need to append it */
  *array_sz = *array_sz + 1;
  *array = (int*) realloc(*array, *array_sz * sizeof(int));
  ADJ_CHKMALLOC(*array);
  (*array)[*array_sz - 1] = value;
  return ADJ_OK;
}

int adj_extend_timestep_data(adj_adjointer* adjointer, int extent)
{
  /* We have an array adjointer->timestep_data, of size adjointer->ntimesteps.
     We want to realloc that to have size extent. We'll also need to zero/initialise
     all the timestep_data's we've just allocated. */
  int i;

  assert(extent > adjointer->ntimesteps);
  adjointer->timestep_data = (adj_timestep_data*) realloc(adjointer->timestep_data, extent * sizeof(adj_timestep_data));
  ADJ_CHKMALLOC(adjointer->timestep_data);
  for (i = adjointer->ntimesteps; i < extent; i++)
  {
    adjointer->timestep_data[i].start_equation = -1;
    adjointer->timestep_data[i].start_time = ADJ_UNSET;
    adjointer->timestep_data[i].end_time = ADJ_UNSET;
    adjointer->timestep_data[i].functional_data_start = NULL;
    adjointer->timestep_data[i].functional_data_end = NULL;
  }
  adjointer->ntimesteps = extent;

  return ADJ_OK;
}

int adj_minval(int* array, int array_sz)
{
  int i;
  int minval;

  assert(array_sz > 0);

  minval = array[0];
  for (i = 1; i < array_sz; i++)
  {
    if (array[i] < minval)
    {
      minval = array[i];
    }
  }

  return minval;
}

/* Provides the number of timesteps that need this variable for the computation of the specified functional */
int adj_variable_get_ndepending_timesteps(adj_adjointer* adjointer, adj_variable variable, char* functional, int* ntimesteps)
{
  int ierr;
  adj_variable_data* data_ptr;
  int k;

  ierr = adj_find_variable_data(&(adjointer->varhash), &variable, &data_ptr);
  if (ierr != ADJ_OK) return ierr;

  *ntimesteps = 0;
  for (k = 0; k < data_ptr->ndepending_timesteps; k++)
  {
    adj_functional_data* func_ptr;
    func_ptr = adjointer->timestep_data[data_ptr->depending_timesteps[k]].functional_data_start;
    while (func_ptr != NULL)
    {
      if (strncmp(func_ptr->name, functional, ADJ_NAME_LEN) == 0)
      {
        /* OK, we've found the functional data for this timestep. Now, we loop through the dependencies of this
           functional at this timestep, to see if we need to add it */
        int l;
        for (l = 0; l < func_ptr->ndepends; l++)
        {
          if (adj_variable_equal(&variable, &func_ptr->dependencies[l], 1))
          {
            *ntimesteps = *ntimesteps + 1;
            break;
          }
        }
        break;
      }
      else
      {
        func_ptr = func_ptr->next;
      }
    }
  }

  return ADJ_OK;
}

/* Provides the timesteps that need this variable for the computation of the specified functional */
int adj_variable_get_depending_timestep(adj_adjointer* adjointer, adj_variable variable, char* functional, int i, int* timestep)
{
  int ierr;
  adj_variable_data* data_ptr;
  int k;
  int ntimesteps;

  ierr = adj_find_variable_data(&(adjointer->varhash), &variable, &data_ptr);
  if (ierr != ADJ_OK) return ierr;

  ntimesteps = 0;
  for (k = 0; k < data_ptr->ndepending_timesteps; k++)
  {
    adj_functional_data* func_ptr;
    func_ptr = adjointer->timestep_data[data_ptr->depending_timesteps[k]].functional_data_start;
    while (func_ptr != NULL)
    {
      if (strncmp(func_ptr->name, functional, ADJ_NAME_LEN) == 0)
      {
        /* OK, we've found the functional data for this timestep. Now, we loop through the dependencies of this
           functional at this timestep, to see if we need to add it */
        int l;
        for (l = 0; l < func_ptr->ndepends; l++)
        {
          if (adj_variable_equal(&variable, &func_ptr->dependencies[l], 1))
          {
            if (i == ntimesteps)
            {
              *timestep = data_ptr->depending_timesteps[k];
              return ADJ_OK;
            }
            else
            {
              ntimesteps = ntimesteps + 1;
              break;
            }
          }
        }
        break;
      }
      else
      {
        func_ptr = func_ptr->next;
      }
    }
  }

  return ADJ_ERR_INVALID_INPUTS;
}
