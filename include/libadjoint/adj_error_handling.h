#ifndef ADJ_ERROR_HANDLING_H
#define ADJ_ERROR_HANDLING_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ADJ_ERROR_MSG_BUF 1024

char adj_error_msg[ADJ_ERROR_MSG_BUF];

#define ADJ_ERR_OK 0
#define ADJ_ERR_INVALID_INPUTS 1
#define ADJ_ERR_HASH_FAILED 2
#define ADJ_ERR_NEED_CALLBACK 3
#define ADJ_ERR_NEED_VALUE 4

char adj_error_codes[3][ADJ_ERROR_MSG_BUF];

#define adj_chkierr(ierr) adj_chkierr_private(ierr, __FILE__, __LINE__)

void adj_init_error_codes(void);
void adj_chkierr_private(int ierr, char* file, int line);
#endif
