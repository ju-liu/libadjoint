#include "libadjoint/adj_error_handling.h"
#include "libadjoint/adj_constants.h"

static int error_check = ADJ_FALSE;

void adj_chkierr_private(int ierr, char* file, int line)
{
  char adj_error_codes[10][ADJ_ERROR_MSG_BUF] = {"ADJ_OK", "ADJ_ERR_INVALID_INPUTS", "ADJ_ERR_HASH_FAILED",
                                                "ADJ_ERR_NEED_CALLBACK", "ADJ_ERR_NEED_VALUE", "ADJ_ERR_NOT_IMPLEMENTED",
                                                "ADJ_ERR_DICT_FAILED", "ADJ_ERR_TOLERANCE_EXCEEDED", "ADJ_ERR_MALLOC_FAILED", 
                                                "ADJ_ERR_REVOLVE_ERROR"};
  char adj_warn_codes[4][ADJ_ERROR_MSG_BUF] = {"ADJ_WARN_ALREADY_RECORDED", "ADJ_WARN_COMPARISON_FAILED", "ADJ_WARN_UNINITIALISED_VALUE",
                                               "ADJ_WARN_NOT_IMPLEMENTED"};

  if (ierr > 0)
  {
    fprintf(stderr, "Error: file %s:%d\n", file, line);
    fprintf(stderr, "Error: got error code %s\n", adj_error_codes[ierr]);
    fprintf(stderr, "Error: %s\n", strlen(adj_error_msg) == 0 ? "(no error message)" : adj_error_msg);
    exit(ierr);
  }
  else if (ierr < 0)
  {
    fprintf(stderr, "Warning: file %s:%d\n", file, line);
    fprintf(stderr, "Warning: got error code %s\n", adj_warn_codes[abs(ierr)-1]);
    fprintf(stderr, "Warning: %s\n", strlen(adj_error_msg) == 0 ? "(no error message)" : adj_error_msg);
  }
}

int adj_chkierr_auto_private(int ierr, char* file, int line)
{
  if (ierr != ADJ_OK && error_check == ADJ_TRUE) adj_chkierr_private(ierr, file, line);
  return ierr;
}

int adj_set_error_checking(int check)
{
  if (check != ADJ_TRUE && check != ADJ_FALSE)
  {
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "check must be either ADJ_TRUE or ADJ_FALSE.");
    return adj_chkierr_auto(ADJ_ERR_INVALID_INPUTS);
  }

  error_check = check;
  return ADJ_OK;
}

