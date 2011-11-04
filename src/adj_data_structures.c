#include "libadjoint/adj_data_structures.h"
#include "libadjoint/adj_error_handling.h"

int adj_create_variable(char* name, int timestep, int iteration, int auxiliary, adj_variable* var)
{
  size_t slen;

  /* zero any struct padding bytes; we'll be hashing this later */
  memset(var, 0, sizeof(adj_variable));

  slen = strlen(name);
  if (slen > ADJ_NAME_LEN)
  {
    strncpy(adj_error_msg, "Name variable too long; recompile with bigger ADJ_NAME_LEN.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (timestep < 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Cannot have a negative timestep %d.", timestep);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (iteration < 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Cannot have a negative iteration %d.", iteration);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  strncpy(var->name, name, ADJ_NAME_LEN);
  var->timestep = timestep;
  var->iteration = iteration;
  var->auxiliary = auxiliary;
  var->type = ADJ_FORWARD;

  return ADJ_OK;
}

int adj_variable_get_name(adj_variable var, char** name)
{
  *name = var.name;
  return ADJ_OK;
}

int adj_variable_get_timestep(adj_variable var, int* timestep)
{
  *timestep = var.timestep;
  return ADJ_OK;
}

int adj_variable_get_iteration(adj_variable var, int* iteration)
{
  *iteration = var.iteration;
  return ADJ_OK;
}

int adj_create_nonlinear_block(char* name, int ndepends, adj_variable* depends, void* context, adj_scalar coefficient, adj_nonlinear_block* nblock)
{
  size_t slen;

  /* zero any struct padding bytes */
  memset(nblock, 0, sizeof(adj_nonlinear_block));

  slen = strlen(name);
  if (slen > ADJ_NAME_LEN)
  {
    strncpy(adj_error_msg, "Name variable too long; recompile with bigger ADJ_NAME_LEN.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (ndepends <= 0)
  {
    strncpy(adj_error_msg, "For it to be nonlinear, it needs at least one dependency.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  strncpy(nblock->name, name, ADJ_NAME_LEN);
  nblock->coefficient = coefficient;
  nblock->context = context;
  nblock->ndepends = ndepends;
  nblock->depends = malloc(ndepends * sizeof(adj_variable));
  ADJ_CHKMALLOC(nblock->depends);
  memcpy(nblock->depends, depends, ndepends * sizeof(adj_variable));
  nblock->test_deriv_hermitian = ADJ_FALSE;
  nblock->number_of_tests = 0;
  nblock->tolerance = (adj_scalar) 0.0;
  nblock->test_derivative = ADJ_FALSE;
  nblock->number_of_rounds = 0;
  return ADJ_OK;
}

int adj_destroy_nonlinear_block(adj_nonlinear_block* nblock)
{
  free(nblock->depends);
  return ADJ_OK;
}

int adj_nonlinear_block_set_coefficient(adj_nonlinear_block* nblock, adj_scalar coefficient)
{
  nblock->coefficient = coefficient;
  return ADJ_OK;
}

int adj_create_block(char* name, adj_nonlinear_block* nblock, void* context, adj_scalar coefficient, adj_block* block)
{
  size_t slen;

  /* zero any struct padding bytes */
  memset(block, 0, sizeof(adj_block));

  slen = strlen(name);
  if (slen > ADJ_NAME_LEN)
  {
    strncpy(adj_error_msg, "Name variable too long; recompile with bigger ADJ_NAME_LEN.", ADJ_ERROR_MSG_BUF);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  strncpy(block->name, name, ADJ_NAME_LEN);

  if (nblock == NULL)
  {
    block->has_nonlinear_block = 0;
  }
  else
  {
    block->has_nonlinear_block = 1;
    block->nonlinear_block = *nblock;
  }

  block->context = context;
  block->hermitian = 0;
  block->coefficient = coefficient;
  block->test_hermitian = ADJ_FALSE;
  block->number_of_tests = 0;
  block->tolerance = (adj_scalar) 0.0;

  return ADJ_OK;
}

int adj_destroy_block(adj_block* block)
{
  /* Dummy to fool the compiler into not printing a warning */
  (void) block;
  return ADJ_OK;
}

int adj_block_set_coefficient(adj_block* block, adj_scalar coefficient)
{
  block->coefficient = coefficient;
  return ADJ_OK;
}

int adj_block_set_hermitian(adj_block* block, int hermitian)
{
  if (hermitian != ADJ_TRUE && hermitian != ADJ_FALSE)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The hermitian argument should either be ADJ_TRUE or ADJ_FALSE.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  block->hermitian = hermitian;
  return ADJ_OK;
}

int adj_block_set_test_hermitian(adj_block* block, int test_hermitian, int number_of_tests, adj_scalar tolerance) 
{
  if (test_hermitian != ADJ_TRUE && test_hermitian != ADJ_FALSE)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The test_hermitian argument should either be ADJ_TRUE or ADJ_FALSE.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (number_of_tests <= 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The number_of_tests argument must be positive.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (tolerance < (adj_scalar) 0.0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The tolerance argument must be nonnegative.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  block->test_hermitian = test_hermitian;
  block->number_of_tests = number_of_tests;
  block->tolerance = tolerance;
  return ADJ_OK;
}

int adj_nonlinear_block_set_test_hermitian(adj_nonlinear_block* nblock, int test_deriv_hermitian, int number_of_tests, adj_scalar tolerance) 
{
  if (test_deriv_hermitian != ADJ_TRUE && test_deriv_hermitian != ADJ_FALSE)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The test_deriv_hermitian argument should either be ADJ_TRUE or ADJ_FALSE.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (number_of_tests <= 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The number_of_tests argument must be positive.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (tolerance < (adj_scalar) 0.0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The tolerance argument must be nonnegative.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  nblock->test_deriv_hermitian = test_deriv_hermitian;
  nblock->number_of_tests = number_of_tests;
  nblock->tolerance = tolerance;
  return ADJ_OK;
}

int adj_nonlinear_block_set_test_derivative(adj_nonlinear_block* nblock, int test_derivative, int number_of_rounds) 
{
  if (test_derivative != ADJ_TRUE && test_derivative != ADJ_FALSE)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The test_derivative argument should either be ADJ_TRUE or ADJ_FALSE.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  if (number_of_rounds <= 1)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "The number_of_rounds argument must be at least 2.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  nblock->test_derivative = test_derivative;
  nblock->number_of_rounds = number_of_rounds;
  return ADJ_OK;
}

int adj_variable_equal(adj_variable* var1, adj_variable* var2, int nvars)
{
  return memcmp(var1, var2, nvars * sizeof(adj_variable)) == 0 ? 1 : 0;
}

int adj_variable_str(adj_variable var, char* name, size_t namelen)
{
  char buf[255];
  memset(buf, 0, 255*sizeof(char));
  memset(name, 0, namelen*sizeof(char));

  switch (var.type)
  {
  case (ADJ_FORWARD):
    snprintf(buf, 255, ":Forward%s", var.auxiliary ? ":Auxiliary" : "");
    break;
  case (ADJ_ADJOINT):
    snprintf(buf, 255, ":Adjoint[%s]", var.functional);
    break;
  case (ADJ_TLM):
    snprintf(buf, 255, ":Sensitivity[%s]", var.functional);
    break;
  default:
    assert(0);
  }
  snprintf(name, namelen, "%s:%d:%d%s", var.name, var.timestep, var.iteration, buf);
  return ADJ_OK;
}

int adj_create_equation(adj_variable var, int nblocks, adj_block* blocks, adj_variable* targets, adj_equation* equation)
{
  int targets_variable;
  int i;

  /* First, let's check the variable isn't auxiliary.
     Auxiliary means we don't solve an equation for it ... */
  if (var.auxiliary)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Cannot register an equation for an auxiliary variable %s.", buf);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  /* Check we have a sane nblocks */
  if (nblocks < 1)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You need at least one block in an equation.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  /* So we haven't seen this variable before. Let's check that the equation actually references this variable.
     Let's also check that no targets are auxiliary */
  targets_variable = 0;
  for (i = 0; i < nblocks; i++)
  {
    if (adj_variable_equal(&var, &(targets[i]), 1))
      targets_variable = 1;

    if (targets[i].auxiliary)
    {
      char buf[ADJ_NAME_LEN];
      adj_variable_str(var, buf, ADJ_NAME_LEN);
      snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Cannot target an auxiliary variable %s.", buf);
      return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
    }
  }

  if (!targets_variable)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Trying to register an equation for %s, but this equation doesn't target this variable.", buf);
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  /* OK. We've done all the sanity checking we can. Let's build the adj_equation. */
  equation->variable = var;
  equation->nblocks = nblocks;
  equation->blocks = (adj_block*) malloc(nblocks * sizeof(adj_block));
  memcpy(equation->blocks, blocks, nblocks * sizeof(adj_block));
  for (i = 0; i < nblocks; i++)
  {
    if (blocks[i].has_nonlinear_block)
    {
      int ierr;
      ierr = adj_copy_nonlinear_block(blocks[i].nonlinear_block, &equation->blocks[i].nonlinear_block);
      if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    }
  }
  equation->targets = (adj_variable*) malloc(nblocks * sizeof(adj_variable));
  memcpy(equation->targets, targets, nblocks * sizeof(adj_variable));

  equation->nrhsdeps = 0;
  equation->rhsdeps = NULL;
  equation->rhs_context = NULL;

  return ADJ_OK;
}

int adj_equation_set_rhs_dependencies(adj_equation* equation, int nrhsdeps, adj_variable* rhsdeps, void* context)
{
  if (nrhsdeps < 0)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You can't have a negative number of dependencies.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  equation->nrhsdeps = nrhsdeps;
  equation->rhsdeps = (adj_variable*) malloc(nrhsdeps * sizeof(adj_variable));
  memcpy(equation->rhsdeps, rhsdeps, nrhsdeps * sizeof(adj_variable));
  equation->rhs_context = context;
  return ADJ_OK;
}

int adj_destroy_equation(adj_equation* equation)
{
  int i;
  int ierr;

  for (i = 0; i < equation->nblocks; i++)
  {
    if (equation->blocks[i].has_nonlinear_block)
    {
      ierr = adj_destroy_nonlinear_block(&equation->blocks[i].nonlinear_block);
      if (ierr != ADJ_OK) return adj_chkierr_auto(ierr);
    }
  }

  free(equation->blocks); equation->blocks = NULL;
  free(equation->targets); equation->targets = NULL;
  if (equation->nrhsdeps > 0)
  {
    free(equation->rhsdeps); equation->rhsdeps = NULL;
  }

  return ADJ_OK;
}

int adj_variable_set_auxiliary(adj_variable* var, int auxiliary)
{
  var->auxiliary = auxiliary;
  return ADJ_OK;
}

int adj_create_nonlinear_block_derivative(adj_adjointer* adjointer, adj_nonlinear_block nblock, adj_variable fwd, adj_vector contraction, int hermitian, adj_nonlinear_block_derivative* deriv)
{
  if (adjointer->callbacks.vec_duplicate == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "We need the ADJ_VEC_DUPLICATE_CB callback to do nonlinear differentiation.");
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }
  if (adjointer->callbacks.vec_axpy == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "We need the ADJ_VEC_AXPY_CB callback to do nonlinear differentiation.");
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }

  deriv->nonlinear_block = nblock;
  deriv->variable = fwd;
  deriv->hermitian = hermitian;

  /* We make a copy because, when we destroy these, the contractions may or may not be 
     vectors stored for other reasons. (The simplification routine to minimise the amount
     of derivatives we have to compute will create new contractions.) So we make it so that
     the nonlinear_block_derivative ALWAYS owns its contraction, so that we can unambiguously
     decide to deallocate it. */
  adjointer->callbacks.vec_duplicate(contraction, &(deriv->contraction));
  adjointer->callbacks.vec_axpy(&(deriv->contraction), (adj_scalar) 1.0, contraction);

  return ADJ_OK;
}

int adj_destroy_nonlinear_block_derivative(adj_adjointer* adjointer, adj_nonlinear_block_derivative* deriv)
{
  if (adjointer->callbacks.vec_destroy == NULL)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "We need the ADJ_VEC_DESTROY_CB callback to destroy adj_nonlinear_block_derivatives.");
    return adj_chkierr_auto(ADJ_ERR_NEED_CALLBACK);
  }

  adjointer->callbacks.vec_destroy(&(deriv->contraction));
  return ADJ_OK;
}

int adj_copy_nonlinear_block(adj_nonlinear_block src, adj_nonlinear_block* dest)
{
  *dest = src; dest->depends = NULL;

  dest->depends = (adj_variable*) malloc(src.ndepends * sizeof(adj_variable));
  ADJ_CHKMALLOC(dest->depends);

  memcpy(dest->depends, src.depends, src.ndepends * sizeof(adj_variable));

  return ADJ_OK;
}
