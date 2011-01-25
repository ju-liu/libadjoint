#ifndef ADJ_DATA_STRUCTURES_H
#define ADJ_DATA_STRUCTURES_H

#include <string.h>
#include "adj_constants.h"
#include "uthash.h"

typedef struct
{
  char name[ADJ_NAME_LEN];
  int timestep; /* what timestep this variable is associated with */
  int iteration; /* what iteration inside the timestep */
  int type; /* forward, adjoint, or tlm */
  int auxiliary; /* is this a real dependency (a variable that is solved for) or auxiliary */
  int functional; /* which functional or parameter is this associated with (adjoint/tlm variables) */
} adj_variable;

typedef struct
{
  void* ptr; /* a pointer to the user's data */
  int klass; /* a field to be set by the user in case adj_vector masks multiple separate types (scalar, vector, etc.) */
} adj_vector;

typedef struct
{ 
  void* ptr; /* a pointer to the user's data */
  int klass; /* a field to be set by the user in case adj_matrix masks multiple separate types (scalar, vector, etc.) */
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
} adj_block;

typedef struct
{
  adj_variable variable;
  int ntargets;
  adj_block* blocks;
  adj_variable* targets;
} adj_equation;

typedef struct
{
  int storage_type;
  int has_value;

  /* for ADJ_STORAGE_MEMORY */
  adj_vector value;

  /* for ADJ_STORAGE_DISK */
  char* filename;

  /* POD, temporal interpolation, ... */
} adj_storage_data;

typedef struct adj_variable_data
{
  int equation;

  int ntargeting_equations;
  int* targeting_equations;

  int ndepending_equations;
  int* depending_equations;

  int nrhs_equations;
  int* rhs_equations;

  int nadjoint_equations;
  int* adjoint_equations;
  adj_storage_data storage;
  struct adj_variable_data* next;
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
  void (*vec_setvalues)(adj_vector *vec, adj_scalar scalars[]);
  void (*vec_getsize)(adj_vector vec, int *sz);
  void (*vec_divide)(adj_vector numerator, adj_vector denominator, adj_vector *output);

  void (*mat_duplicate)(adj_matrix matin, adj_matrix *matout);
  void (*mat_axpy)(adj_matrix *Y, adj_scalar alpha, adj_matrix X);
  void (*mat_destroy)(adj_matrix *mat);
  void (*mat_getvecs)(adj_matrix mat, adj_vector *left);
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
  int nequations;
  adj_equation* equations;

  adj_variable_data_list vardata;
  adj_variable_hash* varhash;

  int options[ADJ_NO_OPTIONS];

  adj_data_callbacks callbacks;
  adj_op_callback_list nonlinear_colouring_sz_list;
  adj_op_callback_list nonlinear_colouring_list;
  adj_op_callback_list nonlinear_action_list;
  adj_op_callback_list nonlinear_derivative_action_list;
  adj_op_callback_list nonlinear_derivative_assembly_list;
  adj_op_callback_list block_action_list;
  adj_op_callback_list block_assembly_list;
} adj_adjointer;

int adj_create_variable(char* name, int timestep, int iteration, int auxiliary, adj_variable* var);
int adj_variable_get_name(adj_variable var, char** name);
int adj_variable_get_timestep(adj_variable var, int* timestep);
int adj_variable_get_iteration(adj_variable var, int* iteration);
int adj_create_nonlinear_block(char* name, int ndepends, adj_variable* depends, adj_scalar coefficient, void* context, adj_nonlinear_block* nblock);
int adj_destroy_nonlinear_block(adj_nonlinear_block* nblock);
int adj_create_block(char* name, adj_nonlinear_block* nblock, void* context, int hermitian, adj_block* block);
int adj_destroy_block(adj_block* block);
int adj_variable_equal(adj_variable* var1, adj_variable* var2, int nvars);
int adj_variable_str(adj_variable var, char* name, size_t namelen);

#endif
