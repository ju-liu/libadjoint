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
    strncpy(adj_error_msg, "Name variable too long; recompile with bigger ADJ_NAMELEN", ADJ_ERROR_MSG_BUF);
    return ADJ_ERR_INVALID_INPUTS;
  }

  strncpy(var->name, name, ADJ_NAMELEN);
  var->timestep = timestep;
  var->iteration = iteration;
  var->auxiliary = auxiliary;

  return ADJ_ERR_OK;
}
