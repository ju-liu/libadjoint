#include "libadjoint/adj_variable_lookup.h"
#include "libadjoint/adj_error_handling.h"

int adj_add_variable_data(adj_variable_hash** hash, adj_variable* var, adj_variable_data* data)
{
  adj_variable_hash* entry;
  adj_variable_hash* check;

  data->type = var->type;

  HASH_FIND(hh, *hash, var, sizeof(adj_variable), check);

  if (check != NULL)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(*var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Hash failed: %s already has an entry.", buf);
    return adj_chkierr_auto(ADJ_ERR_HASH_FAILED);
  }

  entry = (adj_variable_hash*) malloc(sizeof(adj_variable_hash));
  ADJ_CHKMALLOC(entry);
  memset(entry, 0, sizeof(adj_variable_hash));
  entry->variable = *var;
  entry->data = data;

  HASH_ADD(hh, *hash, variable, sizeof(adj_variable), entry);
  return ADJ_OK;
}

int adj_find_variable_data(adj_variable_hash** hash, adj_variable* var, adj_variable_data** data)
{
  adj_variable_hash* check;

  HASH_FIND(hh, *hash, var, sizeof(adj_variable), check);

  if (check == NULL)
  {
    char buf[ADJ_NAME_LEN];
    adj_variable_str(*var, buf, ADJ_NAME_LEN);
    snprintf(adj_error_msg, ADJ_ERROR_MSG_BUF, "Hash failed: %s has no entry.", buf);
    return adj_chkierr_auto(ADJ_ERR_HASH_FAILED);
  }

  *data = check->data;
  return ADJ_OK;
}

void adj_print_hash(adj_variable_hash** hash)
{
  adj_variable_hash* entry;
  adj_variable_hash* tmp;
  HASH_ITER(hh, *hash, entry, tmp)
  {
    char buf[255];
    adj_variable_str(entry->variable, buf, 255); buf[254] = '\0';
    fprintf(stderr, "%s -> %p\n", buf, (void*)entry->data);
  }
}

int adj_destroy_hash(adj_variable_hash** hash)
{
  adj_variable_hash* entry;
  adj_variable_hash* tmp;

  HASH_ITER(hh, *hash, entry, tmp)
  {
    HASH_DEL(*hash, entry);
    free(entry);
  }
  return ADJ_OK;
}
