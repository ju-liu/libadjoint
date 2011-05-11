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

int adj_test_action_transpose(adj_adjointer* adjointer, adj_block block, adj_vector model_input, adj_vector model_output, int N, adj_scalar tol)
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
  adjointer->callbacks.vec_duplicate(x, &ATy);
  adjointer->callbacks.vec_duplicate(model_output, &y);
  adjointer->callbacks.vec_duplicate(y, &Ax);
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

    if (cabs((double complex) yAx - ATyx) > tol) 
    {
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Transpose verification of block \"%s\" failed: |<y, Ax> - <A^Ty, x>| == %lf (> tolerance of %lf).", block.name, cabs((double complex) yAx - ATyx), (double) tol);
      ierr = ADJ_ERR_TOLERANCE_EXCEEDED;
      break;
    }
  }

  adjointer->callbacks.vec_destroy(&x);
  adjointer->callbacks.vec_destroy(&ATy);
  adjointer->callbacks.vec_destroy(&y);
  adjointer->callbacks.vec_destroy(&Ax);

  return ierr;
}
