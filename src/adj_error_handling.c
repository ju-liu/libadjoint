#include "libadjoint/adj_error_handling.h"

void adj_chkierr_private(int ierr, char* file, int line)
{
  char adj_error_codes[7][ADJ_ERROR_MSG_BUF] = {"ADJ_ERR_OK", "ADJ_ERR_INVALID_INPUTS", "ADJ_ERR_HASH_FAILED",
                                                "ADJ_ERR_NEED_CALLBACK", "ADJ_ERR_NEED_VALUE", "ADJ_ERR_NOT_IMPLEMENTED",
                                                "ADJ_ERR_DICT_FAILED"};
  if (ierr > 0)
  {
    fprintf(stderr, "Error: file %s:%d\n", file, line);
    fprintf(stderr, "Error: got error code %s\n", adj_error_codes[ierr]);
    fprintf(stderr, "Error: %s\n", adj_error_msg);
    exit(ierr);
  }
}
