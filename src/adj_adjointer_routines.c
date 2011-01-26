#include "libadjoint/adj_adjointer_routines.h"

/* int adj_create_adjointer(adj_adjointer* adjointer);
int adj_destroy_adjointer(adj_adjointer* adjointer); */

int adj_register_equation(adj_adjointer* adjointer, adj_equation equation)
{
  adj_variable_data data;
  adj_variable_data* new_data;
  int ierr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;

  /* Let's check we haven't solved for this variable before */
  ierr = adj_find_variable_data(adjointer->varhash, &(equation.variable), &data);
  if (ierr != ADJ_ERR_HASH_FAILED)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(equation.variable, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "We have already registered an equation for variable %s.", buf);
    return ADJ_ERR_INVALID_INPUTS;
  }

  /* OK. We're good to go. */
  /* Let's add it to the hash table. */
  new_data = (adj_variable_data*) malloc(sizeof(adj_variable_data));
  new_data->equation = adjointer->nequations + 1; /* we're about to fill it in, don't worry */
  new_data->next = NULL;
  new_data->storage.has_value = 0;

  /* add to the hash table */
  ierr = adj_add_variable_data(adjointer->varhash, &(equation.variable), new_data);
  if (ierr != ADJ_ERR_OK) return ierr;

  /* and add to the data list */
  if (adjointer->vardata.firstnode == NULL)
  {
    adjointer->vardata.firstnode = new_data;
    adjointer->vardata.lastnode = new_data;
  }
  else
  {
    adjointer->vardata.lastnode->next = new_data;
    adjointer->vardata.lastnode = new_data;
  }

  /* OK. Next create an entry for the adj_equation in the adjointer. */

  /* Check we have enough room, and if not, make some */
  if (adjointer->nequations == adjointer->equations_sz)
  {
    adjointer->equations = (adj_equation*) realloc(adjointer->equations, (adjointer->equations_sz + ADJ_PREALLOC_SIZE) * sizeof(adj_equation));
    adjointer->equations_sz = adjointer->equations_sz + ADJ_PREALLOC_SIZE;
  }

  adjointer->nequations++;
  adjointer->equations[adjointer->nequations] = equation;
  /* now we have copies of the pointer to the arrays of targets, blocks, rhs deps. */
  /* but for consistency, any libadjoint object that the user creates, he must destroy --
     it's simpler that way. */
  /* so we're going to make our own copies, so that the user can destroy his. */
  adjointer->equations[adjointer->nequations].blocks = (adj_block*) malloc(equation.nblocks * sizeof(adj_block));
  memcpy(adjointer->equations[adjointer->nequations].blocks, equation.blocks, equation.nblocks * sizeof(adj_block));
  adjointer->equations[adjointer->nequations].targets = (adj_variable*) malloc(equation.nblocks * sizeof(adj_variable));
  memcpy(adjointer->equations[adjointer->nequations].targets, equation.targets, equation.nblocks * sizeof(adj_variable));
  if (equation.nrhsdeps > 0)
  {
    adjointer->equations[adjointer->nequations].rhsdeps = (adj_variable*) malloc(equation.nrhsdeps * sizeof(adj_variable));
    memcpy(adjointer->equations[adjointer->nequations].rhsdeps, equation.rhsdeps, equation.nrhsdeps * sizeof(adj_variable));
  }

  /* Now find all the entries we need to update in the hash table, and update them */

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
  adj_variable_data data;
  int ierr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;

  ierr = adj_find_variable_data(adjointer->varhash, &var, &data);
  if (ierr != ADJ_ERR_OK && !var.auxiliary) return ierr;
  if (ierr != ADJ_ERR_OK && var.auxiliary)
  {
    /* If the variable is auxiliary, it's alright that this is the first time we've ever seen it */
    adj_variable_data* new_data;
    new_data = (adj_variable_data*) malloc(sizeof(adj_variable_data));
    new_data->equation = -1; /* never set */
    new_data->next = NULL;
    new_data->storage.has_value = 0;

    /* add to the hash table */
    ierr = adj_add_variable_data(adjointer->varhash, &var, new_data);
    if (ierr != ADJ_ERR_OK) return ierr;

    /* and add to the data list */
    if (adjointer->vardata.firstnode == NULL)
    {
      adjointer->vardata.firstnode = new_data;
      adjointer->vardata.lastnode = new_data;
    }
    else
    {
      adjointer->vardata.lastnode->next = new_data;
      adjointer->vardata.lastnode = new_data;
    }

    data = *new_data;
  }

  if (data.storage.has_value)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Variable %s already has a value.", buf);
    return ADJ_ERR_INVALID_INPUTS;
  }

  /* Just in case */
  strncpy(adj_error_msg, "Need data callback.", ADJ_ERROR_MSG_BUF);

  switch (storage.storage_type)
  {
    case ADJ_STORAGE_MEMORY:
      if (adjointer->callbacks.vec_duplicate == NULL) return ADJ_ERR_NEED_CALLBACK;
      if (adjointer->callbacks.vec_axpy == NULL) return ADJ_ERR_NEED_CALLBACK;
      data.storage.storage_type = ADJ_STORAGE_MEMORY;
      adjointer->callbacks.vec_duplicate(storage.value, &(data.storage.value));
      adjointer->callbacks.vec_axpy(&(data.storage.value), (adj_scalar)1.0, storage.value);
      break;
    default:
      strncpy(adj_error_msg, "Storage types other than ADJ_STORAGE_MEMORY are not implemented yet.", ADJ_ERROR_MSG_BUF);
      return ADJ_ERR_NOT_IMPLEMENTED;
  }

  return ADJ_ERR_OK;
}


int adj_register_operator_callback(adj_adjointer* adjointer, int type, char* name, void (*fn)(void))
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
      adjointer->callbacks.vec_divide = (void(*)(adj_vector numerator, adj_vector denominator, adj_vector *output)) fn;
      break;
    case ADJ_VEC_SETVALUES_CB:
      adjointer->callbacks.vec_setvalues = (void(*)(adj_vector *vec, adj_scalar scalars[])) fn;
      break;
    case ADJ_VEC_GETSIZE_CB:
      adjointer->callbacks.vec_getsize = (void(*)(adj_vector vec, int *sz)) fn;
      break;

    case ADJ_MAT_DUPLICATE_CB:
      adjointer->callbacks.mat_duplicate = (void(*)(adj_matrix matin, adj_matrix *matout)) fn;
    case ADJ_MAT_AXPY_CB:
      adjointer->callbacks.mat_axpy = (void(*)(adj_matrix *Y, adj_scalar alpha, adj_matrix X)) fn;
    case ADJ_MAT_DESTROY_CB:
      adjointer->callbacks.mat_destroy = (void(*)(adj_matrix *mat)) fn;
    case ADJ_MAT_GETVECS_CB:
      adjointer->callbacks.mat_getvecs = (void(*)(adj_matrix mat, adj_vector *left)) fn;

   default:
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Unknown data callback type %d.", type);
      return ADJ_ERR_INVALID_INPUTS;
  }

  return ADJ_ERR_OK;
}

int adj_forget_adjoint_equation(adj_adjointer* adjointer, int equation)
{
  adj_variable_data* data;
  int should_we_delete;
  int i;
  int ierr;

  if (adjointer->options[ADJ_ACTIVITY] == ADJ_ACTIVITY_NOTHING) return ADJ_ERR_OK;

  data = adjointer->vardata.firstnode;
  while (data != NULL)
  {
    if (data->storage.has_value && data->nadjoint_equations > 0)
    {
      should_we_delete = 1;
      for (i = 0; i < data->nadjoint_equations; i++)
      {
        if (equation > data->adjoint_equations[i])
        {
          should_we_delete = 0;
          break;
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

int adj_get_variable_value(adj_adjointer* adjointer, adj_variable var, adj_vector* value)
{
  int ierr;
  adj_variable_data data;

  ierr = adj_find_variable_data(adjointer->varhash, &var, &data);
  if (ierr != ADJ_ERR_OK) return ierr;

  if (data.storage.storage_type != ADJ_STORAGE_MEMORY)
  {
    ierr = ADJ_ERR_NOT_IMPLEMENTED;
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Sorry, storage strategies other than ADJ_STORAGE_MEMORY are not implemented yet.");
    return ierr;
  }

  if (!data.storage.has_value)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);

    ierr = ADJ_ERR_NEED_VALUE;
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but don't have one recorded.", buf);
    return ierr;
  }

  *value = data.storage.value;
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

adj_storage_data adj_storage_memory(adj_vector value)
{
  adj_storage_data data;

  data.storage_type = ADJ_STORAGE_MEMORY;
  data.value = value;
  return data;
}
