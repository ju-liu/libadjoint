#ifndef HTML_ENCODE_H
#define HTML_ENCODE_H

#include <stdio.h>
#include <string.h>
#include "adj_constants.h"

#ifndef ADJ_HIDE_FROM_USER
char *encode_html(char *html);
char *replace_all_str(char *str, char *orig, char *rep);
char *replace_str(char *str, char *orig, char *rep);
#endif

#endif
