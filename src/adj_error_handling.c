#include "libadjoint/adj_error_handling.h"

void adj_chkierr_private(int ierr, char* file, int line)
{
  if (ierr > 0)
  {
    fprintf(stderr, "Error: file %s:%d\n", file, line);
    fprintf(stderr, "Error: got error code %s\n", adj_error_codes[ierr]);
    fprintf(stderr, "Error: %s\n", adj_error_msg);
    exit(ierr);
  }
}

void adj_init_error_codes(void)
{
  strncpy(adj_error_codes[ADJ_ERR_OK], "ADJ_ERR_OK", ADJ_ERROR_MSG_BUF);
  strncpy(adj_error_codes[ADJ_ERR_INVALID_INPUTS], "ADJ_ERR_INVALID_INPUTS", ADJ_ERROR_MSG_BUF);
}
