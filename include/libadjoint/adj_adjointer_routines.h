#ifndef ADJ_ADJOINTER_ROUTINES_H
#define ADJ_ADJOINTER_ROUTINES_H

#include <assert.h>
#include <stdio.h>
#include "adj_data_structures.h"
#include "adj_variable_lookup.h"
#include "adj_error_handling.h"

int adj_create_adjointer(adj_adjointer* adjointer);
int adj_destroy_adjointer(adj_adjointer* adjointer);
int adj_set_option(adj_adjointer* adjointer, int option, int choice);
int adj_equation_count(adj_adjointer* adjointer, int* count);
int adj_register_equation(adj_adjointer* adjointer, adj_equation equation);
int adj_record_variable(adj_adjointer* adjointer, adj_variable var, adj_storage_data storage);
int adj_register_operator_callback(adj_adjointer* adjointer, int type, char* name, void (*fn)(void));
int adj_register_data_callback(adj_adjointer* adjointer, int type, void (*fn)(void));
int adj_forget_adjoint_equation(adj_adjointer* adjointer, int equation);

int adj_find_operator_callback(adj_adjointer* adjointer, int type, char* name, void (**fn)(void));
int adj_get_variable_value(adj_adjointer* adjointer, adj_variable var, adj_vector* value);
int adj_forget_variable_value(adj_adjointer* adjointer, adj_variable_data* data);
int adj_destroy_variable_data(adj_adjointer* adjointer, adj_variable_data* data);
int adj_add_new_hash_entry(adj_adjointer* adjointer, adj_variable* var, adj_variable_data** data);

adj_storage_data adj_storage_memory(adj_vector value);
void adj_append_unique(int** array, int* array_sz, int value);
#endif
