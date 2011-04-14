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
  adjointer->callbacks.vec_setvalues = NULL;
  adjointer->callbacks.vec_getsize = NULL;
  adjointer->callbacks.vec_divide = NULL;

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
  adjointer->functional_derivative_list.firstnode = NULL;
  adjointer->functional_derivative_list.lastnode = NULL;
  adjointer->forward_source_callback = NULL;

  for (i = 0; i < ADJ_NO_OPTIONS; i++)
    adjointer->options[i] = 0; /* 0 is the default for all options */

  return ADJ_ERR_OK;
}

int adj_destroy_adjointer(adj_adjointer* adjointer)
{
  int i;
  int ierr;
  adj_variable_data* data_ptr;
  adj_variable_data* data_ptr_tmp;
  adj_op_callback* cb_ptr;
  adj_op_callback* cb_ptr_tmp;
  adj_func_deriv_callback* func_deriv_cb_ptr;
  adj_func_deriv_callback* func_deriv_cb_ptr_tmp;
  adj_functional_data* functional_data_ptr_next = NULL;
  adj_functional_data* functional_data_ptr = NULL;

  for (i = 0; i < adjointer->nequations; i++)
  {
    ierr = adj_destroy_equation(&(adjointer->equations[i]));
    if (ierr != ADJ_ERR_OK) return ierr;
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
    if (ierr != ADJ_ERR_OK) return ierr;
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

  func_deriv_cb_ptr = adjointer->functional_derivative_list.firstnode;
  while(func_deriv_cb_ptr != NULL)
  {
    func_deriv_cb_ptr_tmp = func_deriv_cb_ptr;
    func_deriv_cb_ptr = func_deriv_cb_ptr->next;
    free(func_deriv_cb_ptr_tmp);
  }

  adj_create_adjointer(adjointer);
  return ADJ_ERR_OK;
}

int adj_register_equation(adj_adjointer* adjointer, adj_equation equation)
{
  adj_variable_data* data_ptr;
  int ierr;
  int i;
  int j;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;

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

  /* OK. We're good to go. */

  /* Let's add it to the hash table. */
  /* ierr is from the call to adj_find_variable_data. If it is ADJ_ERR_HASH_FAILED, it means
     we have to add it in to the hash table. */
  if (ierr == ADJ_ERR_HASH_FAILED)
  {
    ierr = adj_add_new_hash_entry(adjointer, &(equation.variable), &data_ptr);
    if (ierr != ADJ_ERR_OK) return ierr;
  }
  data_ptr->equation = adjointer->nequations;
  /* OK. Next create an entry for the adj_equation in the adjointer. */

  /* Check we have enough room, and if not, make some */
  if (adjointer->nequations == adjointer->equations_sz)
  {
    adjointer->equations = (adj_equation*) realloc(adjointer->equations, (adjointer->equations_sz + ADJ_PREALLOC_SIZE) * sizeof(adj_equation));
    adjointer->equations_sz = adjointer->equations_sz + ADJ_PREALLOC_SIZE;
  }

  adjointer->nequations++;
  adjointer->equations[adjointer->nequations - 1] = equation;

  /* Do any necessary recording of timestep indices */
  if (adjointer->ntimesteps < equation.variable.timestep + 1) /* adjointer->ntimesteps should be at least equation.variable.timestep + 1 */
  {
    adj_extend_timestep_data(adjointer, equation.variable.timestep + 1); /* extend the array as necessary */
  }
  if (adjointer->timestep_data[equation.variable.timestep].start_equation == -1) /* -1 is the sentinel value for unset */
  {
    adjointer->timestep_data[equation.variable.timestep].start_equation = adjointer->nequations - 1; /* fill in the start equation */
  }

  /* now we have copies of the pointer to the arrays of targets, blocks, rhs deps. */
  /* but for consistency, any libadjoint object that the user creates, he must destroy --
     it's simpler that way. */
  /* so we're going to make our own copies, so that the user can destroy his. */
  adjointer->equations[adjointer->nequations - 1].blocks = (adj_block*) malloc(equation.nblocks * sizeof(adj_block));
  memcpy(adjointer->equations[adjointer->nequations - 1].blocks, equation.blocks, equation.nblocks * sizeof(adj_block));
  adjointer->equations[adjointer->nequations - 1].targets = (adj_variable*) malloc(equation.nblocks * sizeof(adj_variable));
  memcpy(adjointer->equations[adjointer->nequations - 1].targets, equation.targets, equation.nblocks * sizeof(adj_variable));
  if (equation.nrhsdeps > 0)
  {
    adjointer->equations[adjointer->nequations - 1].rhsdeps = (adj_variable*) malloc(equation.nrhsdeps * sizeof(adj_variable));
    memcpy(adjointer->equations[adjointer->nequations - 1].rhsdeps, equation.rhsdeps, equation.nrhsdeps * sizeof(adj_variable));
  }

  /* Now find all the entries we need to update in the hash table, and update them */
  /* First: targeting equations */

  for (i = 0; i < equation.nblocks; i++)
  {
    ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.targets[i]), &data_ptr);
    if (ierr != ADJ_ERR_OK)
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
      if (ierr != ADJ_ERR_OK) return ierr;
      new_data->equation = -1; /* it doesn't have an equation */
      data_ptr = new_data;
    }
    else
    {
      return ierr;
    }

    /* Now data_ptr points to the data we're storing */
    adj_append_unique(&(data_ptr->rhs_equations), &(data_ptr->nrhs_equations), adjointer->nequations - 1);
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
    if (ierr != ADJ_ERR_OK) return ierr;

    eqn_no = data_ptr->equation;
    assert(eqn_no >= 0);

    for (j = 0; j < equation.nrhsdeps; j++)
    {
      ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.rhsdeps[j]), &data_ptr);
      if (ierr != ADJ_ERR_OK) return ierr;
      adj_append_unique(&(data_ptr->adjoint_equations), &(data_ptr->nadjoint_equations), eqn_no); /* dependency j is necessary for equation i */
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
          if (ierr != ADJ_ERR_OK) return ierr;
          new_data->equation = -1; /* it doesn't have an equation */
          data_ptr = new_data;
        }
        else if (ierr == ADJ_ERR_HASH_FAILED)
        {
          return ierr;
        }
        adj_append_unique(&(data_ptr->depending_equations), &(data_ptr->ndepending_equations), adjointer->nequations - 1);
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
      if (ierr != ADJ_ERR_OK) return ierr;

      for (j = 0; j < equation.blocks[i].nonlinear_block.ndepends; j++)
      {
        int k;
        adj_variable_data* j_data;

        /* j_data ALWAYS refers to the data associated with the j'th dependency, throughout this whole loop */
        ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.blocks[i].nonlinear_block.depends[j]), &j_data);
        if (ierr != ADJ_ERR_OK) return ierr;

        /* One set of dependencies: the (adjoint equation of) (the target of this block) (needs) (this dependency) */
        adj_append_unique(&(j_data->adjoint_equations), &(j_data->nadjoint_equations), block_target_data->equation);

        /* Another set of dependencies: the (adjoint equation of) (the j'th dependency) (needs) (the target of this block) */
        adj_append_unique(&(block_target_data->adjoint_equations), &(block_target_data->nadjoint_equations), j_data->equation);

        /* Now we loop over all the dependencies again and fill in the cross-dependencies */
        for (k = 0; k < equation.blocks[i].nonlinear_block.ndepends; k++)
        {
          adj_variable_data* k_data;

          /* k_data ALWAYS refers to the data associated with the k'th dependency, throughout this whole loop */
          ierr = adj_find_variable_data(&(adjointer->varhash), &(equation.blocks[i].nonlinear_block.depends[k]), &k_data);
          if (ierr != ADJ_ERR_OK) return ierr;

          /* Another set of dependencies: the (adjoint equation of) (the j'th dependency) (needs) (the k'th dependency) */
          adj_append_unique(&(k_data->adjoint_equations), &(k_data->nadjoint_equations), j_data->equation);
        }
      }
    }
  }
  return ADJ_ERR_OK;
}

int adj_set_option(adj_adjointer* adjointer, int option, int choice)
{
  if (option < 0 || option >= ADJ_NO_OPTIONS)
  {
    strncpy(adj_error_msg, "Unknown option.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  adjointer->options[option] = choice;
  return ADJ_ERR_OK;
}

int adj_equation_count(adj_adjointer* adjointer, int* count)
{
  *count = adjointer->nequations;
  return ADJ_ERR_OK;
}

int adj_record_variable(adj_adjointer* adjointer, adj_variable var, adj_storage_data storage)
{
  adj_variable_data* data_ptr;
  int ierr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;

  ierr = adj_find_variable_data(&(adjointer->varhash), &var, &data_ptr);
  if (ierr != ADJ_ERR_OK && ierr != ADJ_ERR_HASH_FAILED) return ierr;

  if (ierr == ADJ_ERR_HASH_FAILED)
  {
    /* it's alright that this is the first time we've ever seen it */
    adj_variable_data* new_data;
    ierr = adj_add_new_hash_entry(adjointer, &var, &new_data);
    if (ierr != ADJ_ERR_OK) return ierr;
    new_data->equation = -1; /* it doesn't have an equation */
    data_ptr = new_data;
  }

  assert(data_ptr != NULL);

  if (data_ptr->storage.has_value)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Variable %s already has a value.", buf);
    return ADJ_WARN_ALREADY_RECORDED;
  }

  /* Just in case */
  strncpy(adj_error_msg, "Need data callback.", ADJ_ERROR_MSG_BUF);

  switch (storage.storage_type)
  {
    case ADJ_STORAGE_MEMORY_COPY:
      if (adjointer->callbacks.vec_duplicate == NULL) return ADJ_ERR_NEED_CALLBACK;
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

  return ADJ_ERR_OK;
}


int adj_register_operator_callback(adj_adjointer* adjointer, int type, char* name, void (*fn)(void))
{
  adj_op_callback_list* cb_list_ptr;
  adj_op_callback* cb_ptr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;

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
      return ADJ_ERR_OK;
    }
    cb_ptr = cb_ptr->next;
  }

  /* If we got here, that means that we didn't find it. Tack it on to the end of the list. */
  cb_ptr = (adj_op_callback*) malloc(sizeof(adj_op_callback));
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

  return ADJ_ERR_OK;
}

int adj_register_data_callback(adj_adjointer* adjointer, int type, void (*fn)(void))
{
  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;

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
    case ADJ_VEC_SETVALUES_CB:
      adjointer->callbacks.vec_setvalues = (void(*)(adj_vector *vec, adj_scalar scalars[])) fn;
      break;
    case ADJ_VEC_GETSIZE_CB:
      adjointer->callbacks.vec_getsize = (void(*)(adj_vector vec, int *sz)) fn;
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

  return ADJ_ERR_OK;
}

int adj_register_functional_derivative_callback(adj_adjointer* adjointer, char* name, void (*fn)(adj_adjointer* adjointer, adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output))
{
  adj_func_deriv_callback_list* cb_list_ptr;
  adj_func_deriv_callback* cb_ptr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;

  cb_list_ptr = &(adjointer->functional_derivative_list);

  /* First, we look for an existing callback data structure that might already exist, to replace the function */
  cb_ptr = cb_list_ptr->firstnode;
  while (cb_ptr != NULL)
  {
    if (strncmp(cb_ptr->name, name, ADJ_NAME_LEN) == 0)
    {
      cb_ptr->callback = (void (*)(void* adjointer, adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output)) fn;
      return ADJ_ERR_OK;
    }
    cb_ptr = cb_ptr->next;
  }

  /* If we got here, that means that we didn't find it. Tack it on to the end of the list. */
  cb_ptr = (adj_func_deriv_callback*) malloc(sizeof(adj_func_deriv_callback));
  strncpy(cb_ptr->name, name, ADJ_NAME_LEN);
  cb_ptr->callback = (void (*)(void* adjointer, adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output)) fn;
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

  return ADJ_ERR_OK;
}

int adj_register_forward_source_callback(adj_adjointer* adjointer, void (*fn)(adj_adjointer* adjointer, adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, adj_vector* output))
{
  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;
  adjointer->forward_source_callback = fn;
  return ADJ_ERR_OK;
}

int adj_forget_adjoint_equation(adj_adjointer* adjointer, int equation)
{
  adj_variable_data* data;
  int should_we_delete;
  int i;
  int ierr;
  int min_timestep;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;

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

      /* Also check that it isn't necessary for any timesteps we still have to compute
         functional right-hand-sides for */
      if (data->ndepending_timesteps > 0 && should_we_delete == 1)
      {
        int start_equation;
        min_timestep = adj_minval(data->depending_timesteps, data->ndepending_timesteps);
        ierr = adj_timestep_start_equation(adjointer, min_timestep, &start_equation);
        assert(ierr == ADJ_ERR_OK);

        if (equation > start_equation)
        {
          should_we_delete = 0;
        }
      }

      if (should_we_delete)
      {
        ierr = adj_forget_variable_value(adjointer, data);
        if (ierr != ADJ_ERR_OK) return ierr;
      }
    }

    data = data->next;
  }

  return ADJ_ERR_OK;
}

int adj_find_operator_callback(adj_adjointer* adjointer, int type, char* name, void (**fn)(void))
{
  adj_op_callback_list* cb_list_ptr;
  adj_op_callback* cb_ptr;

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
      return ADJ_ERR_OK;
    }
    cb_ptr = cb_ptr->next;
  }

  snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Could not find callback type %d for operator %s.", type, name);
  return ADJ_ERR_NEED_CALLBACK;
}

int adj_find_functional_derivative_callback(adj_adjointer* adjointer, char* name, void (**fn)(adj_adjointer* adjointer, adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output))
{
  adj_func_deriv_callback_list* cb_list_ptr;
  adj_func_deriv_callback* cb_ptr;

  cb_list_ptr = &(adjointer->functional_derivative_list);

  cb_ptr = cb_list_ptr->firstnode;
  while (cb_ptr != NULL)
  {
    if (strncmp(cb_ptr->name, name, ADJ_NAME_LEN) == 0)
    {
      *fn = (void (*)(adj_adjointer* adjointer, adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output)) cb_ptr->callback;
      return ADJ_ERR_OK;
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
  if (ierr != ADJ_ERR_OK) return ierr;

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
  return ADJ_ERR_OK;
}

int adj_has_variable_value(adj_adjointer* adjointer, adj_variable var)
{
  int ierr;
  adj_variable_data* data_ptr;

  ierr = adj_find_variable_data(&(adjointer->varhash), &var, &data_ptr);
  if (ierr != ADJ_ERR_OK) return ierr;

  if (!data_ptr->storage.has_value)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);

    ierr = ADJ_ERR_NEED_VALUE;
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but don't have one recorded.", buf);
    return ierr;
  }
  return ADJ_ERR_OK;
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
  return ADJ_ERR_OK;
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

  return ADJ_ERR_OK;
}

int adj_storage_memory_copy(adj_vector value, adj_storage_data* data)
{
  data->has_value = 1;
  data->storage_type = ADJ_STORAGE_MEMORY_COPY;
  data->value = value;
  return ADJ_ERR_OK;
}

int adj_storage_memory_incref(adj_vector value, adj_storage_data* data)
{
  data->has_value = 1;
  data->storage_type = ADJ_STORAGE_MEMORY_INCREF;
  data->value = value;
  return ADJ_ERR_OK;
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
  (*data)->ndepending_timesteps = 0;
  (*data)->depending_timesteps = NULL;
  (*data)->nadjoint_equations = 0;
  (*data)->adjoint_equations = NULL;

  /* add to the hash table */
  ierr = adj_add_variable_data(&(adjointer->varhash), var, *data);
  if (ierr != ADJ_ERR_OK) return ierr;

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

  return ADJ_ERR_OK;
}

int adj_timestep_count(adj_adjointer* adjointer, int* count)
{
  *count = adjointer->ntimesteps;
  return ADJ_ERR_OK;
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
  while (ierr == ADJ_ERR_OK);

  if (*count==0)
  {
    strncpy(adj_error_msg, "Error in adj_iteration_count: No iteration found for supplied variable.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }
  return ADJ_ERR_OK;
}

int adj_timestep_start_equation(adj_adjointer* adjointer, int timestep, int* start)
{
  if (timestep < 0 || timestep >= adjointer->ntimesteps)
  {
    strncpy(adj_error_msg, "Invalid timestep supplied to adj_timestep_start.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  *start = adjointer->timestep_data[timestep].start_equation;
  return ADJ_ERR_OK;
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
  return ADJ_ERR_OK;
}

int adj_timestep_set_times(adj_adjointer* adjointer, int timestep, adj_scalar start, adj_scalar end)
{
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
    adj_extend_timestep_data(adjointer, timestep + 1);

  adjointer->timestep_data[timestep].start_time = start;
  adjointer->timestep_data[timestep].end_time = end;

  return ADJ_ERR_OK;
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
    return ADJ_ERR_OK;
}

int adj_timestep_set_functional_dependencies(adj_adjointer* adjointer, int timestep, char* functional, int ndepends, adj_variable* dependencies)
{
  int i;
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
    adj_extend_timestep_data(adjointer, timestep + 1);

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

  if (adjointer->timestep_data[timestep].functional_data_start == NULL)
    adjointer->timestep_data[timestep].functional_data_start = functional_data_ptr;
  if (adjointer->timestep_data[timestep].functional_data_end != NULL) 
    adjointer->timestep_data[timestep].functional_data_end->next = functional_data_ptr;
  adjointer->timestep_data[timestep].functional_data_end = functional_data_ptr;

  functional_data_ptr->next = NULL;
  strncpy(functional_data_ptr->name, functional, ADJ_NAME_LEN);
  functional_data_ptr->ndepends = ndepends;
  functional_data_ptr->dependencies = (adj_variable*) malloc(ndepends * sizeof(adj_variable));
  memcpy(functional_data_ptr->dependencies, dependencies, ndepends * sizeof(adj_variable));

  for (i = 0; i < ndepends; i++)
  {
    int ierr;
    adj_variable_data* data_ptr;

    ierr = adj_find_variable_data(&(adjointer->varhash), &(dependencies[i]), &data_ptr);
    if (ierr == ADJ_ERR_HASH_FAILED)
    {
      ierr = adj_add_new_hash_entry(adjointer, &(dependencies[i]), &data_ptr);
      if (ierr != ADJ_ERR_OK) return ierr;
    }

    /* Record that this variable is necessary for the functional evaluation at this point in time */
    adj_append_unique(&(data_ptr->depending_timesteps), &(data_ptr->ndepending_timesteps), timestep);
  }
  /* We are done */
  return ADJ_ERR_OK;
}

void adj_append_unique(int** array, int* array_sz, int value)
{
  int i;

  for (i = 0; i < *array_sz; i++)
    if ((*array)[i] == value)
      return;

  /* So if we got here, we really do need to append it */
  *array_sz = *array_sz + 1;
  *array = (int*) realloc(*array, *array_sz * sizeof(int));
  (*array)[*array_sz - 1] = value;
  return;
}

void adj_extend_timestep_data(adj_adjointer* adjointer, int extent)
{
  /* We have an array adjointer->timestep_data, of size adjointer->ntimesteps.
     We want to realloc that to have size extent. We'll also need to zero/initialise
     all the timestep_data's we've just allocated. */
  int i;

  assert(extent > adjointer->ntimesteps);
  adjointer->timestep_data = (adj_timestep_data*) realloc(adjointer->timestep_data, extent * sizeof(adj_timestep_data));
  for (i = adjointer->ntimesteps; i < extent; i++)
  {
    adjointer->timestep_data[i].start_equation = -1;
    adjointer->timestep_data[i].start_time = ADJ_UNSET;
    adjointer->timestep_data[i].end_time = ADJ_UNSET;
    adjointer->timestep_data[i].functional_data_start = NULL;
    adjointer->timestep_data[i].functional_data_end = NULL;
  }
  adjointer->ntimesteps = extent;
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

int adj_variable_get_ndepending_timesteps(adj_adjointer* adjointer, adj_variable variable, char* functional, int* ntimesteps)
{
  int ierr;
  adj_variable_data* data_ptr;
  int k;

  ierr = adj_find_variable_data(&(adjointer->varhash), &variable, &data_ptr);
  if (ierr != ADJ_ERR_OK) return ierr;

  *ntimesteps = 0;
  for (k = 0; k < data_ptr->ndepending_timesteps; k++)
  {
    adj_functional_data* func_ptr;
    func_ptr = adjointer->timestep_data[k].functional_data_start;
    while (func_ptr != NULL)
    {
      if (strncmp(func_ptr->name, functional, ADJ_NAME_LEN) == 0)
      {
        *ntimesteps = *ntimesteps + 1;
        break;
      }
      else
      {
        func_ptr = func_ptr->next;
      }
    }
  }

  return ADJ_ERR_OK;
}

int adj_variable_get_depending_timestep(adj_adjointer* adjointer, adj_variable variable, char* functional, int i, int* timestep)
{
  int ierr;
  adj_variable_data* data_ptr;
  int k;
  int ntimesteps;

  ierr = adj_find_variable_data(&(adjointer->varhash), &variable, &data_ptr);
  if (ierr != ADJ_ERR_OK) return ierr;

  ntimesteps = 0;
  for (k = 0; k < data_ptr->ndepending_timesteps; k++)
  {
    adj_functional_data* func_ptr;
    func_ptr = adjointer->timestep_data[k].functional_data_start;
    while (func_ptr != NULL)
    {
      if (strncmp(func_ptr->name, functional, ADJ_NAME_LEN) == 0)
      {
        if (i == ntimesteps)
        {
          *timestep = data_ptr->depending_timesteps[k];
          return ADJ_ERR_OK;
        }
        else
        {
          ntimesteps++;
          break;
        }
      }
      else
      {
        func_ptr = func_ptr->next;
      }
    }
  }

  return ADJ_ERR_INVALID_INPUTS;
}

int adj_adjointer_check_consistency(adj_adjointer* adjointer)
{
  int ierr;
  adj_variable_hash* hash_ptr;

  hash_ptr = adjointer->varhash;
  while(hash_ptr != NULL)
  {
    adj_variable var = hash_ptr->variable;
    adj_variable_data* data_ptr = hash_ptr->data;

    if (var.auxiliary == ADJ_FALSE && data_ptr->equation < 0)
    {
      char buf[ADJ_NAME_LEN];
      adj_variable_str(var, buf, ADJ_NAME_LEN);
      ierr = ADJ_ERR_INVALID_INPUTS;
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Variable %s is not auxiliary, but has no equation set for it.", buf);
      return ierr;
    }

    hash_ptr = hash_ptr->hh.next;
  }

  return ADJ_ERR_OK;
}
