#include "libadjoint/adj_dictionary.h"
#include "libadjoint/adj_error_handling.h"

int adj_dict_init(adj_dictionary* dict)
{
  dict->dict = NULL;
  return ADJ_ERR_OK;
}

int adj_dict_set(adj_dictionary* dict, char* key, char* value)
{
  adj_dictionary_entry* entry;
  adj_dictionary_entry* check;

  HASH_FIND(hh, dict->dict, key, ADJ_DICT_LEN * sizeof(char), check);
  if (check != NULL)
  {
    strncpy(check->value, value, ADJ_DICT_LEN);
    return ADJ_ERR_OK;
  }

  entry = (adj_dictionary_entry*) malloc(sizeof(adj_dictionary_entry));
  memset(entry, 0, sizeof(adj_dictionary_entry));
  strncpy(entry->key, key, ADJ_DICT_LEN * sizeof(char));
  strncpy(entry->value, value, ADJ_DICT_LEN * sizeof(char));

  HASH_ADD(hh, dict->dict, key, ADJ_DICT_LEN * sizeof(char), entry);
  return ADJ_ERR_OK;
}

int adj_dict_find(adj_dictionary* dict, char* key, char** value)
{
  adj_dictionary_entry* check;
  HASH_FIND(hh, dict->dict, key, ADJ_DICT_LEN * sizeof(char), check);

  if (check == NULL) return ADJ_ERR_DICT_FAILED;

  *value = check->value;
  return ADJ_ERR_OK;
}

void adj_dict_print(adj_dictionary* dict)
{
  adj_dictionary_entry* entry;
  adj_dictionary_entry* tmp;
  HASH_ITER(hh, dict->dict, entry, tmp)
  {
    fprintf(stderr, "%s -> %s\n", entry->key, entry->value);
  }
}

int adj_dict_destroy(adj_dictionary* dict)
{
  adj_dictionary_entry* entry;
  adj_dictionary_entry* tmp;
  HASH_ITER(hh, dict->dict, entry, tmp)
  {
    HASH_DEL(dict->dict, entry);
    free(entry);
  }
  return ADJ_ERR_OK;
}