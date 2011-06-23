#include "libadjoint/adj_debug.h"

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

int adj_test_block_action_transpose(adj_adjointer* adjointer, adj_block block, adj_vector model_input, adj_vector model_output, int N, adj_scalar tol)
{
  int i;
  int ierr;
  void (*block_action_func)(int, adj_variable*, adj_vector*, int, adj_scalar, adj_vector, void*, adj_vector*) = NULL;
  adj_vector x, y, Ax, ATy;
  adj_scalar yAx, ATyx;

  if (adjointer->callbacks.vec_set_random == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to test the transpose of an operator, you need the ADJ_VEC_SET_RANDOM_CB callback.");
    return ADJ_ERR_NEED_CALLBACK;
  }
  if (adjointer->callbacks.vec_dot_product == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to test the transpose of an operator, you need the ADJ_VEC_DOT_PRODUCT_CB callback.");
    return ADJ_ERR_NEED_CALLBACK;
  }
  if (adjointer->callbacks.vec_duplicate == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to test the transpose of an operator, you need the ADJ_VEC_DUPLICATE_CB callback.");
    return ADJ_ERR_NEED_CALLBACK;
  }
  if (adjointer->callbacks.vec_destroy == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to test the transpose of an operator, you need the ADJ_VEC_DESTROY_CB callback.");
    return ADJ_ERR_NEED_CALLBACK;
  }

  ierr = adj_find_operator_callback(adjointer, ADJ_BLOCK_ACTION_CB, block.name, (void (**)(void)) &block_action_func);
  if (ierr != ADJ_ERR_OK)
    return ierr;


  adjointer->callbacks.vec_duplicate(model_input, &x);
  adjointer->callbacks.vec_duplicate(model_output, &y);
  block.test_hermitian = ADJ_FALSE;

  for (i = 0; i < N; i++)
  {
    adjointer->callbacks.vec_set_random(&x);
    adjointer->callbacks.vec_set_random(&y);

    ierr = adj_evaluate_block_action(adjointer, block, x, &Ax);
    if (ierr != ADJ_ERR_OK)
      break;

    block.hermitian = !block.hermitian;
    ierr = adj_evaluate_block_action(adjointer, block, y, &ATy);
    if (ierr != ADJ_ERR_OK)
      break;
    block.hermitian = !block.hermitian;

    adjointer->callbacks.vec_dot_product(x, ATy, &ATyx);
    adjointer->callbacks.vec_dot_product(y, Ax, &yAx);
    adjointer->callbacks.vec_destroy(&ATy);
    adjointer->callbacks.vec_destroy(&Ax);

    if (cabs((double complex) yAx - ATyx) > tol) 
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Transpose verification of block \"%s\" failed: |<y, Ax> - <A^Ty, x>| == %e (> tolerance of %e, iteration %d).", block.name, cabs((double complex) yAx - ATyx), (double) tol, i+1);
      ierr = ADJ_ERR_TOLERANCE_EXCEEDED;
      break;
    }
  }

  adjointer->callbacks.vec_destroy(&x);
  adjointer->callbacks.vec_destroy(&y);

  return ierr;
}

int adj_test_nonlinear_derivative_action_transpose(adj_adjointer* adjointer, adj_nonlinear_block_derivative nonlinear_block_derivative, adj_vector model_input, adj_vector model_output, int N, adj_scalar tol)
{
  int i;
  int ierr;
  void (*nonlinear_derivative_action_func)(int nvar, adj_variable* variables, adj_vector* dependencies, adj_variable derivative, adj_vector contraction, int hermitian, adj_vector input, adj_scalar coefficient, void* context, adj_vector* output) = NULL;
  adj_vector x, y, Gx, GTy;
  adj_scalar yGx, GTyx;

  if (adjointer->callbacks.vec_set_random == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to test the transpose of an operator, you need the ADJ_VEC_SET_RANDOM_CB callback.");
    return ADJ_ERR_NEED_CALLBACK;
  }
  if (adjointer->callbacks.vec_dot_product == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to test the transpose of an operator, you need the ADJ_VEC_DOT_PRODUCT_CB callback.");
    return ADJ_ERR_NEED_CALLBACK;
  }
  if (adjointer->callbacks.vec_duplicate == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to test the transpose of an operator, you need the ADJ_VEC_DUPLICATE_CB callback.");
    return ADJ_ERR_NEED_CALLBACK;
  }
  if (adjointer->callbacks.vec_destroy == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "In order to test the transpose of an operator, you need the ADJ_VEC_DESTROY_CB callback.");
    return ADJ_ERR_NEED_CALLBACK;
  }

  ierr = adj_find_operator_callback(adjointer, ADJ_NBLOCK_DERIVATIVE_ACTION_CB, nonlinear_block_derivative.nonlinear_block.name, (void (**)(void)) &nonlinear_derivative_action_func);
  if (ierr != ADJ_ERR_OK)
    return ierr;


  adjointer->callbacks.vec_duplicate(model_input, &x);
  adjointer->callbacks.vec_duplicate(model_output, &y);
  nonlinear_block_derivative.nonlinear_block.test_deriv_hermitian = ADJ_FALSE;

  for (i = 0; i < N; i++)
  {
    adjointer->callbacks.vec_set_random(&x);
    adjointer->callbacks.vec_set_random(&y);
    adjointer->callbacks.vec_duplicate(model_output, &Gx);
    adjointer->callbacks.vec_duplicate(model_input, &GTy);

    ierr = adj_evaluate_nonlinear_derivative_action(adjointer, 1, &nonlinear_block_derivative, x, &Gx);
    if (ierr != ADJ_ERR_OK)
      break;

    nonlinear_block_derivative.hermitian = !nonlinear_block_derivative.hermitian;
    ierr = adj_evaluate_nonlinear_derivative_action(adjointer, 1, &nonlinear_block_derivative, y, &GTy);
    if (ierr != ADJ_ERR_OK)
      break;
    nonlinear_block_derivative.hermitian = !nonlinear_block_derivative.hermitian;

    adjointer->callbacks.vec_dot_product(x, GTy, &GTyx);
    adjointer->callbacks.vec_dot_product(y, Gx, &yGx);
    adjointer->callbacks.vec_destroy(&GTy);
    adjointer->callbacks.vec_destroy(&Gx);

    if (cabs((double complex) yGx - GTyx) > tol) 
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Transpose verification of the derivative of nonlinear block \"%s\" failed: |<y, Gx> - <G^Ty, x>| == %e (> tolerance of %e, iteration %d).", nonlinear_block_derivative.nonlinear_block.name, cabs((double complex) yGx - GTyx), (double) tol, i+1);
      ierr = ADJ_ERR_TOLERANCE_EXCEEDED;
      break;
    }
  }

  adjointer->callbacks.vec_destroy(&x);
  adjointer->callbacks.vec_destroy(&y);

  return ierr;
}

