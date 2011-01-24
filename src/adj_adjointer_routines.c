#include "libadjoint/adj_adjointer_routines.h"

/* int adj_create_adjointer(adj_adjointer* adjointer);
int adj_destroy_adjointer(adj_adjointer* adjointer);
int adj_set_option(adj_adjointer* adjointer, int option, int choice);
int adj_equation_count(adj_adjointer* adjointer, int* count);
int adj_register_equation(adj_adjointer* adjointer, adj_variable var, int nblocks, adj_block* blocks, adj_variable* targets, int nrhsdeps, adj_variable* rhsdeps);
int adj_record_variable(adj_adjointer* adjointer, adj_variable var, adj_vector value);
int adj_record_auxiliary(adj_adjointer* adjointer, adj_variable var, adj_vector value);
int adj_register_operator_callback(adj_adjointer* adjointer, int type, char* name, void (*fn)(void)); */

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
    if (data->has_value && data->nadjoint_equations > 0)
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
    if (strncmp(cb_ptr->name, name, ADJ_NAMELEN) == 0)
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
  if (ierr != ADJ_ERR_OK)
  {
    char buf[255];
    adj_variable_str(var, buf, 255);

    ierr = ADJ_ERR_NEED_VALUE;
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Need a value for %s, but don't have one recorded.", buf);
    return ierr;
  }

  *value = data.value;
  return ADJ_ERR_OK;
}

int adj_forget_variable_value(adj_adjointer* adjointer, adj_variable_data* data)
{
  if (adjointer->callbacks.vec_destroy == NULL)
  {
    strncpy(adj_error_msg, "Need vec_destroy data callback.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_NEED_CALLBACK;
  }
  data->has_value = 0;
  adjointer->callbacks.vec_destroy(&(data->value));
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

  if (data->has_value)
    return adj_forget_variable_value(adjointer, data);

  return ADJ_ERR_OK;
}
