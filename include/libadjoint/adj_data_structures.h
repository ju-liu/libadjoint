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

int adj_create_variable(char* name, int timestep, int iteration, int auxiliary, adj_variable* var)

#endif
