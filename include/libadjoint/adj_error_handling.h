#ifndef ADJ_ERROR_HANDLING_H
#define ADJ_ERROR_HANDLING_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ADJ_ERROR_MSG_BUF 1024

char adj_error_msg[ADJ_ERROR_MSG_BUF];

#define ADJ_OK 0
#define ADJ_ERR_INVALID_INPUTS 1
#define ADJ_ERR_HASH_FAILED 2
#define ADJ_ERR_NEED_CALLBACK 3
#define ADJ_ERR_NEED_VALUE 4
#define ADJ_ERR_NOT_IMPLEMENTED 5
#define ADJ_ERR_DICT_FAILED 6
#define ADJ_ERR_TOLERANCE_EXCEEDED 7
#define ADJ_ERR_MALLOC_FAILED 8
#define ADJ_ERR_REVOLVE_ERROR 9

#define ADJ_WARN_ALREADY_RECORDED -1
#define ADJ_WARN_COMPARISON_FAILED -2
#define ADJ_WARN_UNINITIALISED_VALUE -3
#define ADJ_WARN_NOT_IMPLEMENTED -4

/* if you add a new one, make sure to add it into adj_error_codes in src/adj_error_handling.c */

#define adj_chkierr(ierr) adj_chkierr_private(ierr, __FILE__, __LINE__)
void adj_chkierr_private(int ierr, char* file, int line);

#define adj_chkierr_auto(ierr) adj_chkierr_auto_private(ierr, __FILE__, __LINE__)
int adj_chkierr_auto_private(int ierr, char* file, int line);

int adj_set_error_checking(int check);

#define ADJ_CHKMALLOC(x) \
  if ((void*) (x) == NULL) {\
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Memory allocation failed.");\
    return ADJ_ERR_MALLOC_FAILED;\
  }
#endif
