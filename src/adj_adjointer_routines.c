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

  adjointer->callbacks.mat_duplicate = NULL;
  adjointer->callbacks.mat_axpy = NULL;
  adjointer->callbacks.mat_destroy = NULL;

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

  adjointer->functional_data_start = NULL;
  adjointer->functional_data_end = NULL;

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

  functional_data_ptr = adjointer->functional_data_start;
  while (functional_data_ptr != NULL)
  {
    functional_data_ptr_next = functional_data_ptr->next;
    if (functional_data_ptr->dependencies != NULL) free(functional_data_ptr->dependencies);
    free(functional_data_ptr);
    functional_data_ptr = functional_data_ptr_next;
  }

  if (adjointer->timestep_data != NULL) free(adjointer->timestep_data);

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

int adj_register_equation(adj_adjointer* adjointer, adj_equation equation)
{
  adj_variable_data* data_ptr;
  int ierr;
  int i;
  int j;

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

  if (!data_ptr->storage.has_value) /* If we don't have a value recorded, any compare or overwrite flags can be ignored */
  {
    return adj_record_variable_core(adjointer, data_ptr, storage);
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
      ierr = adj_record_variable_compare(adjointer, data_ptr, var, storage);
      compare_ierr = ierr;
      if (storage.overwrite)
      {
        int record_ierr;

        ierr = adj_forget_variable_value(adjointer, data_ptr);
        if (ierr != ADJ_OK) return ierr;

        record_ierr = adj_record_variable_core(adjointer, data_ptr, storage); /* Overwrite the result anyway */
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
        return adj_record_variable_core(adjointer, data_ptr, storage);
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

int adj_record_variable_core(adj_adjointer* adjointer, adj_variable_data* data_ptr, adj_storage_data storage)
{
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

  switch (storage.storage_type)
  {
    case ADJ_STORAGE_MEMORY_COPY:
      if (adjointer->callbacks.vec_axpy == NULL) return ADJ_ERR_NEED_CALLBACK;
      data_ptr->storage.storage_type = ADJ_STORAGE_MEMORY_COPY;
      data_ptr->storage.has_value = storage.has_value;
      adjointer->callbacks.vec_duplicate(storage.value, &(data_ptr->storage.value));
      adjointer->callbacks.vec_axpy(&(data_ptr->storage.value), (adj_scalar)1.0, storage.value);
      break;
    case ADJ_STORAGE_MEMORY_INCREF:
      data_ptr->storage = storage;
      break;
    default:
      strncpy(adj_error_msg, "Storage types other than ADJ_STORAGE_MEMORY_COPY and ADJ_STORAGE_MEMORY_INCREF  are not implemented yet.", ADJ_ERROR_MSG_BUF);
      return ADJ_ERR_NOT_IMPLEMENTED;
  }

  return ADJ_OK;
}

int adj_record_variable_compare(adj_adjointer* adjointer, adj_variable_data* data_ptr, adj_variable var, adj_storage_data storage)
{
  if (data_ptr->storage.has_value && storage.compare)
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
    if (storage.storage_type != ADJ_STORAGE_MEMORY_COPY && storage.storage_type != ADJ_STORAGE_MEMORY_INCREF)
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

    case ADJ_MAT_DUPLICATE_CB:
      adjointer->callbacks.mat_duplicate = (void(*)(adj_matrix matin, adj_matrix *matout)) fn;
      break;
    case ADJ_MAT_AXPY_CB:
      adjointer->callbacks.mat_axpy = (void(*)(adj_matrix *Y, adj_scalar alpha, adj_matrix X)) fn;
      break;
    case ADJ_MAT_DESTROY_CB:
      adjointer->callbacks.mat_destroy = (void(*)(adj_matrix *mat)) fn;
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

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_OK;

  if (equation >= adjointer->nequations)
  {
    strncpy(adj_error_msg, "No such equation.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  data = adjointer->vardata.firstnode;
  while (data != NULL)
  {
    if (data->storage.has_value)
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

      /* Also check that it isn't necessary for any remaining functional right-hand-sides */
      if (should_we_delete == 1)
      {
        adj_functional_data* functional_data_ptr = NULL;
        adj_variable variable = adjointer->equations[data->equation].variable;

        functional_data_ptr = adjointer->functional_data_start;
        while (functional_data_ptr != NULL)
        { 
          for (i = 0; i < functional_data_ptr->ndepends; i++)
          {
            if (adj_variable_equal(&(functional_data_ptr->dependencies[i]), &variable, 1))
              should_we_delete = 0;
          }
    
          functional_data_ptr = functional_data_ptr->next;
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

int adj_find_functional_callback(adj_adjointer* adjointer, char* name, void (**fn)(adj_adjointer* adjointer, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output))
{
  adj_func_callback_list* cb_list_ptr;
  adj_func_callback* cb_ptr;

  cb_list_ptr = &(adjointer->functional_list);

  cb_ptr = cb_list_ptr->firstnode;
  while (cb_ptr != NULL)
  {
    if (strncmp(cb_ptr->name, name, ADJ_NAME_LEN) == 0)
    {
      *fn = (void (*)(adj_adjointer* adjointer, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output)) cb_ptr->callback;
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

  if (!data_ptr->storage.has_value)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);

    ierr = ADJ_ERR_NEED_VALUE;
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but don't have one recorded.", buf);
    return ierr;
  }

  if ((data_ptr->storage.storage_type != ADJ_STORAGE_MEMORY_COPY) && (data_ptr->storage.storage_type != ADJ_STORAGE_MEMORY_INCREF))
  {
    ierr = ADJ_ERR_NOT_IMPLEMENTED;
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Sorry, storage strategies other than ADJ_STORAGE_MEMORY_COPY and ADJ_STORAGE_MEMORY_INCREF are not implemented yet.");
    return ierr;
  }

  *value = data_ptr->storage.value;
  return ADJ_OK;
}

int adj_has_variable_value(adj_adjointer* adjointer, adj_variable var)
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

  if (!data_ptr->storage.has_value)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);

    ierr = ADJ_ERR_NEED_VALUE;
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but don't have one recorded.", buf);
    return ierr;
  }
  return ADJ_OK;
}

int adj_forget_variable_value(adj_adjointer* adjointer, adj_variable_data* data)
{
  if (adjointer->callbacks.vec_destroy == NULL)
  {
    strncpy(adj_error_msg, "Need vec_destroy data callback.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_NEED_CALLBACK;
  }

  assert(data->storage.has_value);

  data->storage.has_value = 0;
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

  if (data->storage.has_value)
    return adj_forget_variable_value(adjointer, data);

  return ADJ_OK;
}

int adj_storage_memory_copy(adj_vector value, adj_storage_data* data)
{
  data->has_value = ADJ_TRUE;
  data->storage_type = ADJ_STORAGE_MEMORY_COPY;
  data->value = value;
  data->compare = ADJ_FALSE;
  data->comparison_tolerance = (adj_scalar)0.0;
  data->overwrite = ADJ_FALSE;
  return ADJ_OK;
}

int adj_storage_memory_incref(adj_vector value, adj_storage_data* data)
{
  data->has_value = ADJ_TRUE;
  data->storage_type = ADJ_STORAGE_MEMORY_INCREF;
  data->value = value;
  data->compare = ADJ_FALSE;
  data->comparison_tolerance = (adj_scalar)0.0;
  data->overwrite = ADJ_FALSE;
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
  (*data)->storage.has_value = 0;
  (*data)->storage.storage_type = ADJ_UNSET;
  (*data)->ntargeting_equations = 0;
  (*data)->targeting_equations = NULL;
  (*data)->ndepending_equations = 0;
  (*data)->depending_equations = NULL;
  (*data)->nrhs_equations = 0;
  (*data)->rhs_equations = NULL;
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

int adj_set_functional_dependencies(adj_adjointer* adjointer, char* functional, int ndepends, adj_variable* dependencies)
{
  int i;

  adj_functional_data* functional_data_ptr = NULL;

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

  /* Make sure that the dependencies for this functional have not been set before */
  functional_data_ptr = adjointer->functional_data_start;
  while (functional_data_ptr != NULL)
  {
    if (strncmp(functional_data_ptr->name, functional, ADJ_NAME_LEN) == 0)
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The dependencies for functional %s have been already set.", functional);
      return ADJ_ERR_INVALID_INPUTS;
    }
    functional_data_ptr = functional_data_ptr->next;
  }

  /* append the functional information to the adjointer */  
  functional_data_ptr = (adj_functional_data*) malloc(sizeof(adj_functional_data));
  ADJ_CHKMALLOC(functional_data_ptr);

  if (adjointer->functional_data_start == NULL)
    adjointer->functional_data_start = functional_data_ptr;
  if (adjointer->functional_data_end != NULL) 
    adjointer->functional_data_end->next = functional_data_ptr;
  adjointer->functional_data_end = functional_data_ptr;

  functional_data_ptr->next = NULL;
  strncpy(functional_data_ptr->name, functional, ADJ_NAME_LEN);
  functional_data_ptr->ndepends = ndepends;
  functional_data_ptr->dependencies = (adj_variable*) malloc(ndepends * sizeof(adj_variable));
  ADJ_CHKMALLOC(functional_data_ptr->dependencies);
  memcpy(functional_data_ptr->dependencies, dependencies, ndepends * sizeof(adj_variable));
  functional_data_ptr->functional_derivative_data_start = NULL;
  functional_data_ptr->functional_derivative_data_end = NULL;

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

  }
  /* We are done */
  return ADJ_OK;
}

int adj_set_functional_derivative_dependencies(adj_adjointer* adjointer, char* functional, adj_variable derivative, int ndepends, adj_variable* dependencies)
{
  int i;

  adj_functional_data* functional_data_ptr = NULL;
  adj_functional_derivative_data* functional_derivative_data_ptr = NULL;

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

  if (derivative.type != ADJ_FORWARD)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Argument derivative must be forward variable.");
    return ADJ_ERR_INVALID_INPUTS;
  }

  /* Find the functional data associated with this functional */
  functional_data_ptr = adjointer->functional_data_start;
  while (functional_data_ptr != NULL)
  {
    if (strncmp(functional_data_ptr->name, functional, ADJ_NAME_LEN) == 0)
      break;
    functional_data_ptr = functional_data_ptr->next;
  }

  if (functional_data_ptr == NULL) 
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The functional dependencies for %s have to be defined first before the derivative dependencies can be set.", functional);
    return ADJ_ERR_INVALID_INPUTS;
  }

  /* Make sure that the dependencies for this functional have not been set before */
  functional_derivative_data_ptr = functional_data_ptr->functional_derivative_data_start;
  while (functional_derivative_data_ptr != NULL)
  {
    if (adj_variable_equal(&(functional_derivative_data_ptr->derivative), &derivative, 1))
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The dependencies for derivative of functional %s with respect to variable %s have been already set.", functional, derivative.name);
      return ADJ_ERR_INVALID_INPUTS;
    }
    functional_derivative_data_ptr = functional_derivative_data_ptr->next;
  }

  /* append the functional derivative dependency information to the list */  
  functional_derivative_data_ptr = (adj_functional_derivative_data*) malloc(sizeof(adj_functional_derivative_data));
  ADJ_CHKMALLOC(functional_derivative_data_ptr);

  if (functional_data_ptr->functional_derivative_data_start == NULL)
    functional_data_ptr->functional_derivative_data_start = functional_derivative_data_ptr;
  if (functional_data_ptr->functional_derivative_data_end != NULL) 
    functional_data_ptr->functional_derivative_data_end->next = functional_derivative_data_ptr;
  functional_data_ptr->functional_derivative_data_end = functional_derivative_data_ptr;

  functional_derivative_data_ptr->next = NULL;
  functional_derivative_data_ptr->ndepends = ndepends;
  functional_derivative_data_ptr->dependencies = (adj_variable*) malloc(ndepends * sizeof(adj_variable));
  ADJ_CHKMALLOC(functional_derivative_data_ptr->dependencies);
  memcpy(functional_derivative_data_ptr->dependencies, dependencies, ndepends * sizeof(adj_variable));

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


