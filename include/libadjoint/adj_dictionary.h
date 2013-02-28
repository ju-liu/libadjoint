#ifndef ADJ_DICTIONARY_H
#define ADJ_DICTIONARY_H

#include "adj_data_structures.h"

#ifdef __cplusplus
extern "C" {
#endif

int adj_dict_init(adj_dictionary* dict);
int adj_dict_set(adj_dictionary* dict, char* key, char* value);
int adj_dict_find(adj_dictionary* dict, char* key, char** value);
void adj_dict_print(adj_dictionary* dict);
int adj_dict_destroy(adj_dictionary* dict);

#ifdef __cplusplus
}
#endif

#endif
