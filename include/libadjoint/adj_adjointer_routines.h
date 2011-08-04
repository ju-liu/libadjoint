#ifndef ADJ_ADJOINTER_ROUTINES_H
#define ADJ_ADJOINTER_ROUTINES_H

#include <assert.h>
#include <stdio.h>
#include "adj_data_structures.h"
#include "adj_variable_lookup.h"
#include "adj_error_handling.h"

int adj_create_adjointer(adj_adjointer* adjointer);
int adj_destroy_adjointer(adj_adjointer* adjointer);
int adj_deactivate_adjointer(adj_adjointer* adjointer);
int adj_equation_count(adj_adjointer* adjointer, int* count);
int adj_register_equation(adj_adjointer* adjointer, adj_equation equation);
int adj_record_variable(adj_adjointer* adjointer, adj_variable var, adj_storage_data storage);
int adj_register_operator_callback(adj_adjointer* adjointer, int type, char* name, void (*fn)(void));
int adj_register_data_callback(adj_adjointer* adjointer, int type, void (*fn)(void));
int adj_register_functional_callback(adj_adjointer* adjointer, char* name, void (*fn)(adj_adjointer* adjointer, int timestep, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output));
int adj_register_functional_derivative_callback(adj_adjointer* adjointer, char* name, void (*fn)(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output));
int adj_register_forward_source_callback(adj_adjointer* adjointer, void (*fn)(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, void* context, adj_vector* output, int* has_output));
int adj_forget_adjoint_equation(adj_adjointer* adjointer, int equation);

int adj_timestep_count(adj_adjointer* adjointer, int* count);
int adj_iteration_count(adj_adjointer* adjointer, adj_variable variable, int* count);
int adj_timestep_start_equation(adj_adjointer* adjointer, int timestep, int* start);
int adj_timestep_end_equation(adj_adjointer* adjointer, int timestep, int* end);
int adj_timestep_set_times(adj_adjointer* adjointer, int timestep, adj_scalar start, adj_scalar end);
int adj_timestep_get_times(adj_adjointer* adjointer, int timestep, adj_scalar* start, adj_scalar* end);
int adj_set_functional_dependencies(adj_adjointer* adjointer, char* functional, int ndepends, adj_variable* dependencies);

int adj_storage_memory_copy(adj_vector value, adj_storage_data* data);
int adj_storage_memory_incref(adj_vector value, adj_storage_data* data);
int adj_storage_set_compare(adj_storage_data* data, int compare, adj_scalar comparison_tolerance);
int adj_storage_set_overwrite(adj_storage_data* data, int overwrite);

#ifndef ADJ_HIDE_FROM_USER
int adj_set_option(adj_adjointer* adjointer, int option, int choice);

int adj_find_operator_callback(adj_adjointer* adjointer, int type, char* name, void (**fn)(void));
int adj_find_functional_callback(adj_adjointer* adjointer, char* name, void (**fn)(adj_adjointer* adjointer, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output));
int adj_find_functional_derivative_callback(adj_adjointer* adjointer, char* name, void (**fn)(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output));
int adj_get_variable_value(adj_adjointer* adjointer, adj_variable var, adj_vector* value);
int adj_has_variable_value(adj_adjointer* adjointer, adj_variable var);
int adj_forget_variable_value(adj_adjointer* adjointer, adj_variable_data* data);
int adj_destroy_variable_data(adj_adjointer* adjointer, adj_variable_data* data);
int adj_add_new_hash_entry(adj_adjointer* adjointer, adj_variable* var, adj_variable_data** data);
int adj_record_variable_core(adj_adjointer* adjointer, adj_variable_data* data_ptr, adj_storage_data storage);
int adj_record_variable_compare(adj_adjointer* adjointer, adj_variable_data* data_ptr, adj_variable var, adj_storage_data storage);

int adj_append_unique(int** array, int* array_sz, int value);
int adj_extend_timestep_data(adj_adjointer* adjointer, int extent);
int adj_extend_functional_data(adj_timestep_data* timestep_data, int extent);
int adj_minval(int* array, int array_sz);
#endif
#endif
