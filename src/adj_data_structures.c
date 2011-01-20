#include "libadjoint/adj_data_structures.h"
#include "libadjoint/adj_error_handling.h"

int adj_create_variable(char* name, int timestep, int iteration, int auxiliary, adj_variable* var)
{
  size_t slen;

  slen = strlen(name);
  if (slen > ADJ_NAMELEN)
  {
    adj_error_msg = "Name variable too long; recompile with bigger ADJ_NAMELEN";
    return ADJ_ERR_INVALID_INPUTS;
  }

  return ADJ_ERR_OK;
}
