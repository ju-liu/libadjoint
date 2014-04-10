#include "libadjoint/html_encode.h"

/* Replaces non html characters with their html code */
char *encode_html(char *html)
{
  static char buf[ADJ_NAME_LEN];
  char *replace_all_buf;

  strncpy(buf, html, ADJ_NAME_LEN);
  replace_all_buf = replace_all_str(buf, "&", "&amp;");
  strncpy(buf, replace_all_buf, ADJ_NAME_LEN);
  replace_all_buf = replace_all_str(buf, "<", "&lt;");
  strncpy(buf, replace_all_buf, ADJ_NAME_LEN);
  replace_all_buf = replace_all_str(buf, ">", "&gt;");
  strncpy(buf, replace_all_buf, ADJ_NAME_LEN);
  replace_all_buf = replace_all_str(buf, "\"", "&quot;");
  strncpy(buf, replace_all_buf, ADJ_NAME_LEN);
  replace_all_buf = replace_all_str(buf, "'", "&apos;");
  strncpy(buf, replace_all_buf, ADJ_NAME_LEN);
  buf[ADJ_NAME_LEN-1] = '\0';

  return buf;
}

/* Replaces all occurrences of orig in str with rep */
char *replace_all_str(char *str, char *orig, char *rep)
{
  static char buf[ADJ_NAME_LEN];
  char *replace_buf, *p;

  strncpy(buf, str, ADJ_NAME_LEN);
  while((p = strstr(buf, orig))) /* Loop as long as orig is in buf */
  {
    replace_buf = replace_str(buf, orig, rep);
    strncpy(buf, replace_buf, ADJ_NAME_LEN);
    buf[ADJ_NAME_LEN-1] = '\0';
  }
  return buf;
}

/* Replaces the first occurrence of orig in str with rep */
char *replace_str(char *str, char *orig, char *rep)
{
  static char buf[ADJ_NAME_LEN];
  char *p;

  if(!(p = strstr(str, orig)))  /* Is 'orig' even in 'str'? */
    return buf;

  strncpy(buf, str, p-str); /* Copy characters from 'str' start to 'orig' st */
  buf[p-str] = '\0';

  sprintf(buf+(p-str), "%s%s", rep, p+strlen(orig));

  return buf;
}
