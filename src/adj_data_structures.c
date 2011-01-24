#include "libadjoint/adj_data_structures.h"
#include "libadjoint/adj_error_handling.h"

int adj_create_variable(char* name, int timestep, int iteration, int auxiliary, adj_variable* var)
{
  size_t slen;

  /* zero any struct padding bytes; we'll be hashing this later */
  memset(var, 0, sizeof(adj_variable));

  slen = strlen(name);
  if (slen > ADJ_NAMELEN)
  {
    strncpy(adj_error_msg, "Name variable too long; recompile with bigger ADJ_NAMELEN.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  strncpy(var->name, name, ADJ_NAMELEN);
  var->timestep = timestep;
  var->iteration = iteration;
  var->auxiliary = auxiliary;

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

int adj_create_nonlinear_block(char* name, int ndepends, adj_variable* depends, adj_scalar coefficient, void* context, adj_nonlinear_block* nblock)
{
  size_t slen;

  /* zero any struct padding bytes */
  memset(nblock, 0, sizeof(adj_nonlinear_block));

  slen = strlen(name);
  if (slen > ADJ_NAMELEN)
  {
    strncpy(adj_error_msg, "Name variable too long; recompile with bigger ADJ_NAMELEN.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  if (ndepends <= 0)
  {
    strncpy(adj_error_msg, "For it to be nonlinear, it needs at least one dependency.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  strncpy(nblock->name, name, ADJ_NAMELEN);
  nblock->coefficient = coefficient;
  nblock->context = context;
  nblock->depends = malloc(ndepends * sizeof(adj_variable));
  memcpy(nblock->depends, depends, ndepends * sizeof(adj_variable));
  return ADJ_ERR_OK;
}

int adj_destroy_nonlinear_block(adj_nonlinear_block* nblock)
{
  free(nblock->depends);
  return ADJ_ERR_OK;
}

int adj_create_block(char* name, adj_nonlinear_block* nblock, void* context, int hermitian, adj_block* block)
{
  size_t slen;

  /* zero any struct padding bytes */
  memset(block, 0, sizeof(adj_block));

  slen = strlen(name);
  if (slen > ADJ_NAMELEN)
  {
    strncpy(adj_error_msg, "Name variable too long; recompile with bigger ADJ_NAMELEN.", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  strncpy(block->name, name, ADJ_NAMELEN);

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
  block->hermitian = hermitian;

  return ADJ_ERR_OK;
}

int adj_destroy_block(adj_block* block)
{
  return ADJ_ERR_OK;
}

int adj_variable_equal(adj_variable* var1, adj_variable* var2, int nvars)
{
  return memcmp(var1, var2, nvars * sizeof(adj_variable)) == 0 ? 1 : 0;
}

int adj_variable_str(adj_variable var, char* name, size_t namelen)
{
  char buf[255];
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
  }
  snprintf(name, namelen, "%s:%d:%d%s", var.name, var.timestep, var.iteration, buf);
  return ADJ_ERR_OK;
}
