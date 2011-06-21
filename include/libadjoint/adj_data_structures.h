#ifndef ADJ_DATA_STRUCTURES_H
#define ADJ_DATA_STRUCTURES_H

#include <string.h>
#include <assert.h>
#include "adj_constants.h"
#include "uthash.h"

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
} adj_nonlinear_block;

typedef struct
{
  char name[ADJ_NAME_LEN];
  int has_nonlinear_block;
  adj_nonlinear_block nonlinear_block;
  void* context;
  int hermitian;
  adj_scalar coefficient;
  int test_hermitian;
  int number_of_tests;
  adj_scalar tolerance;
} adj_block;

typedef struct
{
  adj_variable variable;
  int nblocks;
  adj_block* blocks;
  adj_variable* targets;
  int nrhsdeps;
  adj_variable* rhsdeps;
  void* rhs_context;
} adj_equation;

typedef struct
{
  int storage_type;
  int has_value;

  /* Should we compare against something we already have? */
  int compare;
  adj_scalar comparison_tolerance;

  /* Should we overwrite something that's already recorded? */
  int overwrite;

  /* for ADJ_STORAGE_MEMORY */
  adj_vector value;

  /* for ADJ_STORAGE_DISK */
  char* filename;

  /* POD, temporal interpolation, ... */
} adj_storage_data;

typedef struct adj_variable_data
{
  int equation; /* the equation that solves for this variable */

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
  adj_variable_data* firstnode;
  adj_variable_data* lastnode;
} adj_variable_data_list;

typedef struct
{
  void (*vec_duplicate)(adj_vector x, adj_vector *newx);
  void (*vec_axpy)(adj_vector *y, adj_scalar alpha, adj_vector x);
  void (*vec_destroy)(adj_vector *x);
  void (*vec_set_values)(adj_vector *vec, adj_scalar scalars[]);
  void (*vec_get_size)(adj_vector vec, int *sz);
  void (*vec_divide)(adj_vector *numerator, adj_vector denominator);
  void (*vec_get_norm)(adj_vector x, adj_scalar* norm);
  void (*vec_dot_product)(adj_vector x, adj_vector y, adj_scalar* val);
  void (*vec_set_random)(adj_vector* x);

  void (*mat_duplicate)(adj_matrix matin, adj_matrix *matout);
  void (*mat_axpy)(adj_matrix *Y, adj_scalar alpha, adj_matrix X);
  void (*mat_destroy)(adj_matrix *mat);
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
  void (*callback)(void* adjointer, int timestep, int nb_variables, adj_variable* variables, adj_vector* dependencies, char* name, adj_scalar* output);
  struct adj_func_callback* next;
} adj_func_callback;

typedef struct adj_func_deriv_callback
{
  char name[ADJ_NAME_LEN];
  /* we want this to be adj_adjointer* adjointer, but we haven't defined adj_adjointer yet */
  void (*callback)(void* adjointer, adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output);
  struct adj_func_deriv_callback* next;
} adj_func_deriv_callback;

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
  adj_nonlinear_block nonlinear_block;
  adj_variable variable;
  adj_vector contraction;
  int hermitian;
} adj_nonlinear_block_derivative;

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

typedef struct adj_adjointer
{
  int nequations; /* Number of equations we have registered */
  int equations_sz; /* Number of equations we can store without mallocing -- not the same! */
  adj_equation* equations; /* Array of equations we have registered */

  int ntimesteps; /* Number of timesteps we have seen */
  adj_timestep_data* timestep_data; /* Data for each timestep we have seen */

  adj_variable_hash* varhash; /* The hash table for looking up information about variables */
  adj_variable_data_list vardata; /* We also store a linked list so we can walk all our variable data */

  int options[ADJ_NO_OPTIONS]; /* Pretty obvious */

  adj_data_callbacks callbacks; /* Data callbacks */
  adj_op_callback_list nonlinear_colouring_list; /* Operator callbacks */
  adj_op_callback_list nonlinear_action_list;
  adj_op_callback_list nonlinear_derivative_action_list;
  adj_op_callback_list nonlinear_derivative_assembly_list;
  adj_op_callback_list block_action_list;
  adj_op_callback_list block_assembly_list;
  adj_func_callback_list functional_list;
  adj_func_deriv_callback_list functional_derivative_list;
  void (*forward_source_callback)(struct adj_adjointer* adjointer, adj_variable variable, int nb_variables, adj_variable* variables, adj_vector* dependencies, void* context, adj_vector* output, int* has_output);
} adj_adjointer;

int adj_create_variable(char* name, int timestep, int iteration, int auxiliary, adj_variable* var);
int adj_variable_get_name(adj_variable var, char** name);
int adj_variable_get_timestep(adj_variable var, int* timestep);
int adj_variable_get_iteration(adj_variable var, int* iteration);
int adj_variable_set_auxiliary(adj_variable* var, int auxiliary);
int adj_create_nonlinear_block(char* name, int ndepends, adj_variable* depends, void* context, adj_nonlinear_block* nblock);
int adj_destroy_nonlinear_block(adj_nonlinear_block* nblock);
int adj_nonlinear_block_set_coefficient(adj_nonlinear_block* nblock, adj_scalar coefficient);
int adj_create_block(char* name, adj_nonlinear_block* nblock, void* context, adj_block* block);
int adj_destroy_block(adj_block* block);
int adj_block_set_coefficient(adj_block* block, adj_scalar coefficient);
int adj_block_set_hermitian(adj_block* block, int hermitian);
int adj_block_set_test_hermitian(adj_block* block, int test_hermitian, int number_of_tests, adj_scalar tolerance);
int adj_create_equation(adj_variable var, int nblocks, adj_block* blocks, adj_variable* targets, adj_equation* equation);
int adj_equation_set_rhs_dependencies(adj_equation* equation, int nrhsdeps, adj_variable* rhsdeps, void* context);
int adj_destroy_equation(adj_equation* equation);

#ifndef ADJ_HIDE_FROM_USER
int adj_variable_equal(adj_variable* var1, adj_variable* var2, int nvars);
int adj_variable_str(adj_variable var, char* name, size_t namelen);
int adj_create_nonlinear_block_derivative(adj_adjointer* adjointer, adj_nonlinear_block nblock, adj_variable fwd, adj_vector contraction, int hermitian, adj_nonlinear_block_derivative* deriv);
int adj_destroy_nonlinear_block_derivative(adj_adjointer* adjointer, adj_nonlinear_block_derivative* deriv);
#endif

#endif
