#ifndef ADJ_ERROR_HANDLING_H
#define ADJ_ERROR_HANDLING_H

#include <stdio.h>
#include <stdlib.h>

char adj_error_msg[1024];

#define ADJ_ERR_OK 0
#define ADJ_ERR_INVALID_INPUTS 1

char adj_error_codes[2][1024];

#define adj_chkierr(ierr) adj_chkierr_private(ierr, __FILE__, __LINE__)

void adj_init_error_codes(void);
void adj_chkierr_private(int ierr, char* file, int line);
#endif
