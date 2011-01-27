#include "libadjoint/adj_variable_lookup.h"
#include "libadjoint/adj_error_handling.h"

int adj_add_variable_data(adj_variable_hash** hash, adj_variable* var, adj_variable_data* data)
{
  adj_variable_hash* entry;
  adj_variable_hash* check;

  HASH_FIND(hh, *hash, var, sizeof(adj_variable), check);

  if (check != NULL)
    return ADJ_ERR_HASH_FAILED;

  entry = (adj_variable_hash*) malloc(sizeof(adj_variable_hash));
  memset(entry, 0, sizeof(adj_variable_hash));
  entry->variable = *var;
  entry->data = data;

  HASH_ADD(hh, *hash, variable, sizeof(adj_variable), entry);
  return ADJ_ERR_OK;
}

int adj_find_variable_data(adj_variable_hash** hash, adj_variable* var, adj_variable_data** data)
{
  adj_variable_hash* check;

  HASH_FIND(hh, *hash, var, sizeof(adj_variable), check);

  if (check == NULL)
    return ADJ_ERR_HASH_FAILED;

  *data = check->data;
  return ADJ_ERR_OK;
}

void adj_print_hash(adj_variable_hash** hash)
{
  adj_variable_hash* entry;
  adj_variable_hash* tmp;
  HASH_ITER(hh, *hash, entry, tmp)
  {
    char buf[255];
    adj_variable_str(entry->variable, buf, 255); buf[254] = '\0';
    fprintf(stderr, "%s -> %p\n", buf, entry->data);
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
  return ADJ_ERR_OK;
}
