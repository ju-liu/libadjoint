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
    return ADJ_ERR_INVALID_INPUTS;
  }

  strncpy(var->name, name, ADJ_NAME_LEN);
  var->timestep = timestep;
  var->iteration = iteration;
  var->auxiliary = auxiliary;
  var->type = ADJ_FORWARD;

  return ADJ_ERR_OK;
}

int adj_variable_get_name(adj_variable var, char** name)
{
  *name = var.name;
  return ADJ_ERR_OK;
}

int adj_variable_get_timestep(adj_variable var, int* timestep)
{
  *timestep = var.timestep;
  return ADJ_ERR_OK;
}

int adj_variable_get_iteration(adj_variable var, int* iteration)
{
  *iteration = var.iteration;
  return ADJ_ERR_OK;
}

int adj_create_nonlinear_block(char* name, int ndepends, adj_variable* depends, void* context, adj_nonlinear_block* nblock)
{
  size_t slen;

  /* zero any struct padding bytes */
  memset(nblock, 0, sizeof(adj_nonlinear_block));

  slen = strlen(name);
  if (slen > ADJ_NAME_LEN)
  {
    strncpy(adj_error_msg, "Name variable too long; recompile with bigger ADJ_NAME_LEN.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  if (ndepends <= 0)
  {
    strncpy(adj_error_msg, "For it to be nonlinear, it needs at least one dependency.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  strncpy(nblock->name, name, ADJ_NAME_LEN);
  nblock->coefficient = (adj_scalar)1.0;
  nblock->context = context;
  nblock->ndepends = ndepends;
  nblock->depends = malloc(ndepends * sizeof(adj_variable));
  memcpy(nblock->depends, depends, ndepends * sizeof(adj_variable));
  return ADJ_ERR_OK;
}

int adj_destroy_nonlinear_block(adj_nonlinear_block* nblock)
{
  free(nblock->depends);
  return ADJ_ERR_OK;
}

int adj_nonlinear_block_set_coefficient(adj_nonlinear_block* nblock, adj_scalar coefficient)
{
  nblock->coefficient = coefficient;
  return ADJ_ERR_OK;
}

int adj_create_block(char* name, adj_nonlinear_block* nblock, void* context, adj_block* block)
{
  size_t slen;

  /* zero any struct padding bytes */
  memset(block, 0, sizeof(adj_block));

  slen = strlen(name);
  if (slen > ADJ_NAME_LEN)
  {
    strncpy(adj_error_msg, "Name variable too long; recompile with bigger ADJ_NAME_LEN.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
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

  return ADJ_ERR_OK;
}

int adj_destroy_block(adj_block* block)
{
  /* Dummy to fool the compiler into not printing a warning */
  (void) block;
  return ADJ_ERR_OK;
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
    snprintf(buf, 255, ":Adjoint[%d]", var.functional);
    break;
  case (ADJ_TLM):
    snprintf(buf, 255, ":Sensitivity[%d]", var.functional);
    break;
  default:
    assert(0);
  }
  snprintf(name, namelen, "%s:%d:%d%s", var.name, var.timestep, var.iteration, buf);
  return ADJ_ERR_OK;
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
    return ADJ_ERR_INVALID_INPUTS;
  }

  /* Check we have a sane nblocks */
  if (nblocks < 1)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "You need at least one block in an equation.");
    return ADJ_ERR_INVALID_INPUTS;
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
      return ADJ_ERR_INVALID_INPUTS;
    }
  }

  if (!targets_variable)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Trying to register an equation for %s, but this equation doesn't target this variable.", buf);
    return ADJ_ERR_INVALID_INPUTS;
  }

  /* OK. We've done all the sanity checking we can. Let's build the adj_equation. */
  equation->variable = var;
  equation->nblocks = nblocks;
  equation->blocks = (adj_block*) malloc(nblocks * sizeof(adj_block));
  memcpy(equation->blocks, blocks, nblocks * sizeof(adj_block));
  equation->targets = (adj_variable*) malloc(nblocks * sizeof(adj_variable));
  memcpy(equation->targets, targets, nblocks * sizeof(adj_variable));

  equation->nrhsdeps = 0;
  equation->rhsdeps = NULL;

  return ADJ_ERR_OK;
}

int adj_set_rhs_dependencies(adj_equation* equation, int nrhsdeps, adj_variable* rhsdeps)
{
  if (nrhsdeps < 1)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "If you are registering rhs dependencies, you must have at least one dependency.");
    return ADJ_ERR_INVALID_INPUTS;
  }

  equation->nrhsdeps = nrhsdeps;
  equation->rhsdeps = (adj_variable*) malloc(nrhsdeps * sizeof(adj_variable));
  memcpy(equation->rhsdeps, rhsdeps, nrhsdeps * sizeof(adj_variable));
  return ADJ_ERR_OK;
}

int adj_destroy_equation(adj_equation* equation)
{
  free(equation->blocks); equation->blocks = NULL;
  free(equation->targets); equation->targets = NULL;
  if (equation->nrhsdeps > 0)
  {
    free(equation->rhsdeps); equation->rhsdeps = NULL;
  }

  return ADJ_ERR_OK;
}
