#ifndef ADJ_DATA_STRUCTURES_H
#define ADJ_DATA_STRUCTURES_H

#include <string.h>
#include "adj_constants.h"

typedef struct
{
  char name[ADJ_NAMELEN];
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
  char name[ADJ_NAMELEN];
  adj_scalar coefficient;
  int ndepends;
  adj_variable* depends;
  void* context;
} adj_nonlinear_block;

typedef struct
{
  char name[ADJ_NAMELEN];
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
  int has_value;
  adj_vector value;
  struct adj_variable_data* next;
} adj_variable_data;

typedef struct
{
  adj_variable_data* firstnode;
  adj_variable_data* lastnode;
} adj_variable_data_list;

typedef struct
{
  void (*vec_duplicate)(void);
  void (*vec_axpy)(void);
  void (*vec_destroy)(void);
  void (*vec_setvalues)(void);
  void (*vec_divide)(void);

  void (*mat_duplicate)(void);
  void (*mat_axpy)(void);
  void (*mat_destroy)(void);
  void (*mat_getvecs)(void);
} adj_data_callbacks;

typedef struct adj_op_callback
{
  char name[ADJ_NAMELEN];
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
  int nequations;
  adj_equation* equations;

  adj_variable_data_list vardata;
  void* varhash;

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

#endif
