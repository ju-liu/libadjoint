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
  int nadjoint_equations;
  int* adjoint_equations;
  int has_value;
  adj_vector value;
  struct adj_variable_data* next;
} adj_variable_data;

int adj_create_variable(char* name, int timestep, int iteration, int auxiliary, adj_variable* var);

#endif
