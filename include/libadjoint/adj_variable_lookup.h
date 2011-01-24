#ifndef ADJ_VARIABLE_LOOKUP_H
#define ADJ_VARIABLE_LOOKUP_H

#include "adj_data_structures.h"

int adj_add_variable_data(adj_variable_hash* hash, adj_variable* var, adj_variable_data* data);
int adj_find_variable_data(adj_variable_hash* hash, adj_variable* var, adj_variable_data* data);
int adj_destroy_hash(adj_variable_hash* hash);

#endif
