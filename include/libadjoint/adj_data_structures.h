#ifndef ADJ_DATA_STRUCTURES_H
#define ADJ_DATA_STRUCTURES_H

#include <string.h>
#include <assert.h>
#include "adj_constants.h"
#include "uthash.h"
#include "revolve_c.h"

typedef struct
{
  char name[ADJ_NAME_LEN];
  int timestep; /* what timestep this variable is associated with */
  int iteration; /* what iteration inside the timestep */
  int type; /* forward, adjoint, or tlm: ADJ_FORWARD, ADJ_ADJOINT, or ADJ_TLM */
  int auxiliary; /* is this a real dependency (a variable that is solved for) or auxiliary */
  char functional[ADJ_NAME_LEN]; /* which functional or parameter is this associated with (adjoint/tlm variables) */
} adj_variable;

typedef struct
{
  void* ptr; /* a pointer to the user's data */
  int klass; /* a field to be set by the user in case adj_vector masks multiple separate types (scalar, vector, etc.) */
  int flags; /* for any flags the user might like to set */
} adj_vector;

typedef struct
{ 
  void* ptr; /* a pointer to the user's data */
  int klass; /* a field to be set by the user in case adj_matrix masks multiple separate types (scalar, vector, etc.) */
  int flags; /* for any flags the user might like to set */
} adj_matrix;

typedef struct
{
  char name[ADJ_NAME_LEN];
  adj_scalar coefficient;
  int ndepends;
  adj_variable* depends;
  void* context;
  int test_deriv_hermitian; /* if ADJ_TRUE, the hermitian implementation of this nonlinear block's derivatives are tested automatically */
  int number_of_tests;
  adj_scalar tolerance;
  int test_derivative; /* Flags for test_derivative */
  int number_of_rounds;
} adj_nonlinear_block;

typedef struct
{
  char name[ADJ_NAME_LEN];
  int has_nonlinear_block;
  adj_nonlinear_block nonlinear_block;
  void* context;
  int hermitian;
  adj_scalar coefficient;
  int test_hermitian; /* Flags for test_hermitian */
  int number_of_tests;
  adj_scalar tolerance;
} adj_block;

typedef struct
{
  int nblocks;
  adj_block* blocks;
  adj_variable* targets;
} adj_term;

typedef struct
{
  adj_variable variable;
  int nblocks;
  adj_block* blocks;
  adj_variable* targets;
  int nrhsdeps;
  adj_variable* rhsdeps;
  void* rhs_context;
  void (*rhs_callback)(void* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, void* context, adj_vector* output, int* has_output);
  void (*rhs_deriv_action_callback)(void* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, \
                                    adj_variable d_variable, adj_vector contraction, int hermitian, void* context, adj_vector* output, int* has_output);
  void (*rhs_deriv_assembly_callback)(void* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, \
                                    int hermitian, void* context, adj_matrix* output);
  int memory_checkpoint; /* Can we restart the computation from this equation using variables in memory? */
  int disk_checkpoint; /* Can we restart the computation from this equation using variables on disk? */
} adj_equation;

typedef struct
{
  /* Should we compare against something we already have? */
  int compare;
  adj_scalar comparison_tolerance;

  /* Should we overwrite something that's already recorded? */
  int overwrite;

  adj_vector value;
  /* for ADJ_STORAGE_MEMORY */
  int storage_memory_type; /* ADJ_STORAGE_MEMORY_COPY or ADJ_STORAGE_MEMORY_INCREF */
  int storage_memory_has_value;
  int storage_memory_is_checkpoint; /* memory checkpoints are not deleted by adj_forget_forward_equation */

  /* for ADJ_STORAGE_DISK */
  int storage_disk_has_value;
  int storage_disk_is_checkpoint; /* disk checkpoints are not deleted by adj_forget_forward_equation */

  /* POD, temporal interpolation, ... */
} adj_storage_data;

typedef struct adj_variable_data
{
  int equation; /* the equation that solves for this variable. If the data belongs to a adjoint variable, this will be set to -1 */
  int type; /* is it ADJ_FORWARD, ADJ_ADJOINT or ADJ_TLM? */

  int ntargeting_equations; /* any equations that target this variable */
  int* targeting_equations;

  int ndepending_equations; /* any equations whose operators depend on this variable */
  int* depending_equations;

  int nrhs_equations; /* any equations whose right-hand sides depend on this variable */
  int* rhs_equations;

  int ndepending_timesteps; /* any timesteps that need this variable for the computation of a functional */
  int* depending_timesteps;

  int nadjoint_equations; /* computed: the adjoint equations that need this variable */
  int* adjoint_equations;

  adj_storage_data storage; /* its storage record */
  struct adj_variable_data* next; /* a pointer to the next one, so we can walk the list */
} adj_variable_data;

typedef struct
{
  void (*vec_duplicate)(adj_vector x, adj_vector *newx);
  void (*vec_axpy)(adj_vector *y, adj_scalar alpha, adj_vector x);
  void (*vec_destroy)(adj_vector *x);
  void (*vec_set_values)(adj_vector *vec, adj_scalar scalars[]);
  void (*vec_get_values)(adj_vector vec, adj_scalar *scalars[]);
  void (*vec_get_size)(adj_vector vec, int *sz);
  void (*vec_divide)(adj_vector *numerator, adj_vector denominator);
  void (*vec_get_norm)(adj_vector x, adj_scalar* norm);
  void (*vec_dot_product)(adj_vector x, adj_vector y, adj_scalar* val);
  void (*vec_set_random)(adj_vector* x);
  void (*vec_write)(adj_variable var, adj_vector x);
  void (*vec_read)(adj_variable var, adj_vector* x);
  void (*vec_delete)(adj_variable var);

  void (*mat_duplicate)(adj_matrix matin, adj_matrix *matout);
  void (*mat_axpy)(adj_matrix *Y, adj_scalar alpha, adj_matrix X);
  void (*mat_destroy)(adj_matrix *mat);
  void (*mat_action)(adj_matrix mat, adj_vector x, adj_vector* y);

  void (*solve)(adj_variable var, adj_matrix mat, adj_vector rhs, adj_vector *soln);
} adj_data_callbacks;

typedef struct adj_op_callback
{
  char name[ADJ_NAME_LEN];
  void (*callback)(void);
  struct adj_op_callback* next;
} adj_op_callback;

typedef struct
{
  adj_op_callback* firstnode;
  adj_op_callback* lastnode;
} adj_op_callback_list;

typedef struct adj_func_callback
{
  char name[ADJ_NAME_LEN];
  /* we want this to be adj_adjointer* adjointer, but we haven't defined adj_adjointer yet */
  void (*callback)(void* adjointer, int timestep, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output);
  struct adj_func_callback* next;
} adj_func_callback;

typedef struct adj_func_deriv_callback
{
  char name[ADJ_NAME_LEN];
  /* we want this to be adj_adjointer* adjointer, but we haven't defined adj_adjointer yet */
  void (*callback)(void* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, char* functional, adj_vector* output);
  struct adj_func_deriv_callback* next;
} adj_func_deriv_callback;

typedef struct adj_func_second_deriv_callback
{
  char name[ADJ_NAME_LEN];
  /* we want this to be adj_adjointer* adjointer, but we haven't defined adj_adjointer yet */
  void (*callback)(void* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, adj_vector contraction, char* functional, adj_vector* output);
  struct adj_func_second_deriv_callback* next;
} adj_func_second_deriv_callback;

typedef struct adj_parameter_source_callback
{
  char name[ADJ_NAME_LEN];
  /* we want this to be adj_adjointer* adjointer, but we haven't defined adj_adjointer yet */
  void (*callback)(void* adjointer, int equation, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, char* parameter, adj_vector* output, int* has_output);
  struct adj_parameter_source_callback* next;
} adj_parameter_source_callback;

typedef struct
{
  adj_func_callback* firstnode;
  adj_func_callback* lastnode;
} adj_func_callback_list;

typedef struct
{
  adj_func_deriv_callback* firstnode;
  adj_func_deriv_callback* lastnode;
} adj_func_deriv_callback_list;

typedef struct
{
  adj_func_second_deriv_callback* firstnode;
  adj_func_second_deriv_callback* lastnode;
} adj_func_second_deriv_callback_list;

typedef struct
{
  adj_parameter_source_callback* firstnode;
  adj_parameter_source_callback* lastnode;
} adj_parameter_source_callback_list;

typedef struct
{
  adj_nonlinear_block nonlinear_block; /* nonlinear operator to differentiate */
  adj_variable variable; /* variable to differentiate with respect to */
  adj_vector contraction; /* contraction vector to perform rank reduction */
  int hermitian;
} adj_nonlinear_block_derivative;

typedef struct
{
  adj_nonlinear_block nonlinear_block; /* nonlinear operator to differentiate */
  adj_variable inner_variable; /* the variable for the first derivative */
  adj_vector inner_contraction; /* the contraction for the first derivative */
  adj_variable outer_variable; /* the variable for the second derivative */
  adj_vector outer_contraction; /* the contraction for the second derivative */
  int hermitian;
  adj_vector block_action; /* the variable the second derivative acts on */
} adj_nonlinear_block_second_derivative; /* this structure is needed in the second-order adjoint equation */

typedef struct
{
  adj_variable variable;
  adj_variable_data* data;
  adj_hash_handle hh;
} adj_variable_hash;

typedef struct
{
  char key[ADJ_DICT_LEN];
  char value[ADJ_DICT_LEN];
  adj_hash_handle hh;
} adj_dictionary_entry;

typedef struct
{
  adj_dictionary_entry* dict;
} adj_dictionary;

typedef struct adj_functional_data
{
  char name[ADJ_NAME_LEN];
  int ndepends;
  adj_variable* dependencies;
  struct adj_functional_data* next; /* a pointer to the next one, so we can walk the list */
} adj_functional_data;

typedef struct
{
  int start_equation;

  adj_scalar start_time;
  adj_scalar end_time;

  adj_functional_data* functional_data_start;
  adj_functional_data* functional_data_end;
} adj_timestep_data;

typedef struct
{
  CRevolve revolve; /* The C wrapper of the Revolve object */
  int snaps; /* The total number of available checkpoint slots */
  int snaps_in_ram; /* The number of available RAM checkpoint slots */
  int steps; /* The total number of timesteps in the simulation */
  CACTION current_action; /* The revolve action which is currently being executed */
  int current_timestep; /* The timestep which is currently being executed in the model. */
                        /* This should be thought of as "Libadjoint has all the data to solve */
                        /* for the first forward variable with timestep current_timestep" */
  int verbose; /* A flag that prints revolve specific information to the screen if set to ADJ_TRUE */
  int overwrite; /* A flag indicating if a replay should be performed even if that variable is already recorded. */
                 /* The new value is compared with the existing one in order to check if the revolve replay produces the same solution than the original forward system */
  adj_scalar comparison_tolerance; /* The comparison tolerance in case that overwrite is ADJ_TRUE */
} adj_revolve_data;

typedef struct adj_adjointer
{
  adj_equation* equations; /* Array of equations we have registered */
  int nequations; /* Number of equations we have registered */
  int equations_sz; /* Number of equations we can store without mallocing -- not the same! */

  int ntimesteps; /* Number of timesteps we have seen */
  adj_timestep_data* timestep_data; /* Data for each timestep we have seen */
  
  adj_revolve_data revolve_data; /* A data struct for revolve related information */

  adj_variable_hash* varhash; /* The hash table for looking up information about variables */

  int options[ADJ_NO_OPTIONS]; /* Pretty obvious */

  adj_data_callbacks callbacks; /* Data callbacks */
  adj_op_callback_list nonlinear_colouring_list; /* Operator callbacks */
  adj_op_callback_list nonlinear_action_list;
  adj_op_callback_list nonlinear_derivative_action_list;
  adj_op_callback_list nonlinear_derivative_assembly_list;
  adj_op_callback_list block_action_list;
  adj_op_callback_list block_assembly_list;
  adj_op_callback_list nonlinear_second_derivative_action_list;
  adj_func_callback_list functional_list;
  adj_func_deriv_callback_list functional_derivative_list;
  adj_func_second_deriv_callback_list functional_second_derivative_list;
  adj_parameter_source_callback_list parameter_source_list;

  int finished; /* Is the annotation finished? */
} adj_adjointer;

int adj_create_variable(char* name, int timestep, int iteration, int auxiliary, adj_variable* var);
int adj_variable_get_name(adj_variable var, char** name);
int adj_variable_get_timestep(adj_variable var, int* timestep);
int adj_variable_get_iteration(adj_variable var, int* iteration);
int adj_variable_get_type(adj_variable var, int* type);
int adj_variable_set_auxiliary(adj_variable* var, int auxiliary);
int adj_variable_str(adj_variable var, char* name, size_t namelen);
int adj_create_nonlinear_block(char* name, int ndepends, adj_variable* depends, void* context, adj_scalar coefficient, adj_nonlinear_block* nblock);
int adj_destroy_nonlinear_block(adj_nonlinear_block* nblock);
int adj_nonlinear_block_set_coefficient(adj_nonlinear_block* nblock, adj_scalar coefficient);
int adj_create_block(char* name, adj_nonlinear_block* nblock, void* context, adj_scalar coefficient, adj_block* block);
int adj_destroy_block(adj_block* block);
int adj_block_set_coefficient(adj_block* block, adj_scalar coefficient);
int adj_block_set_hermitian(adj_block* block, int hermitian);
int adj_block_set_test_hermitian(adj_block* block, int test_hermitian, int number_of_tests, adj_scalar tolerance);
int adj_nonlinear_block_set_test_hermitian(adj_nonlinear_block* nblock, int test_deriv_hermitian, int number_of_tests, adj_scalar tolerance); 
int adj_nonlinear_block_set_test_derivative(adj_nonlinear_block* nblock, int test_derivative, int number_of_rounds);
int adj_create_equation(adj_variable var, int nblocks, adj_block* blocks, adj_variable* targets, adj_equation* equation);
int adj_equation_set_rhs_dependencies(adj_equation* equation, int nrhsdeps, adj_variable* rhsdeps, void* context);
int adj_destroy_equation(adj_equation* equation);
int adj_create_term(int nblocks, adj_block* blocks, adj_variable* targets, adj_term* term);
int adj_add_terms(adj_term termA, adj_term termB, adj_term* termC);
int adj_destroy_term(adj_term* term);
int adj_add_term_to_equation(adj_term term, adj_equation* equation);
int adj_equation_set_rhs_callback(adj_equation* equation, void (*fn)(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, void* context, adj_vector* output, int* has_output));
int adj_equation_set_rhs_derivative_action_callback(adj_equation* equation, void (*fn)(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, \
                                    adj_variable d_variable, adj_vector contraction, int hermitian, void* context, adj_vector* output, int* has_output));
int adj_equation_set_rhs_derivative_assembly_callback(adj_equation* equation, void (*fn)(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, \
                                    int hermitian, void* context, adj_matrix* output));
int adj_variable_equal(adj_variable* var1, adj_variable* var2, int nvars);


#ifndef ADJ_HIDE_FROM_USER
int adj_create_nonlinear_block_derivative(adj_adjointer* adjointer, adj_nonlinear_block nblock, adj_scalar block_coefficient, adj_variable fwd, adj_vector contraction, int hermitian, adj_nonlinear_block_derivative* deriv);
int adj_destroy_nonlinear_block_derivative(adj_adjointer* adjointer, adj_nonlinear_block_derivative* deriv);
int adj_create_nonlinear_block_second_derivative(adj_adjointer* adjointer, adj_nonlinear_block nblock, adj_scalar block_coefficient, adj_variable inner_var, adj_vector inner_contraction, adj_variable outer_var, adj_vector outer_contraction, int hermitian, adj_vector action, adj_nonlinear_block_second_derivative* deriv);
int adj_destroy_nonlinear_block_second_derivative(adj_adjointer* adjointer, adj_nonlinear_block_second_derivative* deriv);
int adj_copy_nonlinear_block(adj_nonlinear_block src, adj_nonlinear_block* dest);
int adj_equation_rhs_nonlinear_index(adj_equation eqn);
#endif

#endif
