#ifndef ADJ_ADJOINTER_ROUTINES_H
#define ADJ_ADJOINTER_ROUTINES_H

#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include "adj_data_structures.h"
#include "adj_variable_lookup.h"
#include "adj_error_handling.h"
#include "revolve_c.h"

int adj_create_adjointer(adj_adjointer* adjointer);
int adj_destroy_adjointer(adj_adjointer* adjointer);
int adj_deactivate_adjointer(adj_adjointer* adjointer);
int adj_set_checkpoint_strategy(adj_adjointer* adjointer, int strategy);
int adj_set_revolve_options(adj_adjointer* adjointer, int steps, int snaps_on_disk, int snaps_in_ram, int verbose);
int adj_set_revolve_debug_options(adj_adjointer* adjointer, int overwrite, adj_scalar comparison_tolerance);
int adj_equation_count(adj_adjointer* adjointer, int* count);
int adj_register_equation(adj_adjointer* adjointer, adj_equation equation, int* checkpoint_storage);
int adj_record_variable(adj_adjointer* adjointer, adj_variable var, adj_storage_data storage);
int adj_register_operator_callback(adj_adjointer* adjointer, int type, char* name, void (*fn)(void));
int adj_register_data_callback(adj_adjointer* adjointer, int type, void (*fn)(void));
int adj_register_functional_callback(adj_adjointer* adjointer, char* name, void (*fn)(adj_adjointer* adjointer, int timestep, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output));
int adj_register_functional_derivative_callback(adj_adjointer* adjointer, char* name, void (*fn)(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output));
int adj_register_parameter_source_callback(adj_adjointer* adjointer, char* name, void (*fn)(adj_adjointer* adjointer, int equation, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output, int* has_output));

int adj_forget_adjoint_equation(adj_adjointer* adjointer, int equation);
int adj_forget_forward_equation(adj_adjointer* adjointer, int equation);
int adj_forget_tlm_equation(adj_adjointer* adjointer, int equation);

int adj_forget_adjoint_values(adj_adjointer* adjointer, int equation);
int adj_forget_tlm_values(adj_adjointer* adjointer, int equation);

int adj_timestep_count(adj_adjointer* adjointer, int* count);
int adj_iteration_count(adj_adjointer* adjointer, adj_variable variable, int* count);
int adj_timestep_start_equation(adj_adjointer* adjointer, int timestep, int* start);
int adj_timestep_end_equation(adj_adjointer* adjointer, int timestep, int* end);
int adj_timestep_set_times(adj_adjointer* adjointer, int timestep, adj_scalar start, adj_scalar end);
int adj_timestep_get_times(adj_adjointer* adjointer, int timestep, adj_scalar* start, adj_scalar* end);
int adj_timestep_set_functional_dependencies(adj_adjointer* adjointer, int timestep, char* functional, int ndepends, adj_variable* dependencies);

int adj_storage_memory_copy(adj_vector value, adj_storage_data* data);
int adj_storage_memory_incref(adj_vector value, adj_storage_data* data);
int adj_storage_disk(adj_vector value, adj_storage_data* data);
int adj_storage_set_compare(adj_storage_data* data, int compare, adj_scalar comparison_tolerance);
int adj_storage_set_overwrite(adj_storage_data* data, int overwrite);
int adj_storage_set_checkpoint(adj_storage_data* data, int checkpoint);

int adj_variable_known(adj_adjointer* adjointer, adj_variable var, int* known);
int adj_get_variable_value(adj_adjointer* adjointer, adj_variable var, adj_vector* value);

int adj_set_finished(adj_adjointer* adjointer, int  finished);
int adj_get_finished(adj_adjointer* adjointer, int* finished);

int adj_get_forward_variable(adj_adjointer* adjointer, int i, adj_variable* fwd_var);

#ifndef ADJ_HIDE_FROM_USER
int adj_set_storage_memory_copy(adj_adjointer* adjointer, adj_variable* var);
int adj_set_storage_memory_incref(adj_adjointer* adjointer, adj_variable* var);
int adj_set_option(adj_adjointer* adjointer, int option, int choice);
int adj_variable_get_ndepending_timesteps(adj_adjointer* adjointer, adj_variable variable, char* functional, int* ntimesteps);
int adj_variable_get_depending_timestep(adj_adjointer* adjointer, adj_variable variable, char* functional, int k, int* timestep);
int adj_forget_forward_equation_until(adj_adjointer* adjointer, int equation, int last_equation);
int adj_get_checkpoint_strategy(adj_adjointer* adjointer, int* strategy);

int adj_find_operator_callback(adj_adjointer* adjointer, int type, char* name, void (**fn)(void));
int adj_find_functional_callback(adj_adjointer* adjointer, char* name, void (**fn)(adj_adjointer* adjointer, int timestep, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output));
int adj_find_functional_derivative_callback(adj_adjointer* adjointer, char* functional, void (**fn)(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output));
int adj_find_parameter_source_callback(adj_adjointer* adjointer, char* parameter, void (**fn)(adj_adjointer* adjointer, int equation, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output, int* has_output));
int adj_has_variable_value(adj_adjointer* adjointer, adj_variable var);
int adj_has_variable_value_memory(adj_adjointer* adjointer, adj_variable var);
int adj_has_variable_value_disk(adj_adjointer* adjointer, adj_variable var);
int adj_is_variable_memory_checkpoint(adj_adjointer* adjointer, adj_variable var);
int adj_is_variable_disk_checkpoint(adj_adjointer* adjointer, adj_variable var);
int adj_forget_variable_value(adj_adjointer* adjointer, adj_variable var, adj_variable_data* data);
int adj_forget_variable_value_from_memory(adj_adjointer* adjointer, adj_variable_data* data);
int adj_forget_variable_value_from_disk(adj_adjointer* adjointer, adj_variable var, adj_variable_data* data);
int adj_destroy_variable_data(adj_adjointer* adjointer, adj_variable var, adj_variable_data* data);
int adj_add_new_hash_entry(adj_adjointer* adjointer, adj_variable* var, adj_variable_data** data);
int adj_record_variable_core_disk(adj_adjointer* adjointer, adj_variable var, adj_variable_data* data_ptr, adj_storage_data storage);
int adj_record_variable_core_memory(adj_adjointer* adjointer, adj_variable_data* data_ptr, adj_storage_data storage);
int adj_record_variable_compare(adj_adjointer* adjointer, adj_variable_data* data_ptr, adj_variable var, adj_storage_data storage);

int adj_append_unique(int** array, int* array_sz, int value);
int adj_extend_timestep_data(adj_adjointer* adjointer, int extent);
int adj_extend_functional_data(adj_timestep_data* timestep_data, int extent);
int adj_minval(int* array, int array_sz);
int adj_get_revolve_checkpoint_storage(adj_adjointer* adjointer, adj_equation equation, int* checkpoint_storage); 
int adj_initialise_revolve(adj_adjointer* adjointer);

int adj_checkpoint_equation(adj_adjointer* adjointer, int eqn_number, int checkpoint_strategy);
int adj_checkpoint_variable(adj_adjointer* adjointer, adj_variable var, int checkpoint_strategy);
#endif
#endif
