#ifndef ADJ_PETSC_DATA_STRUCTURES_H
#define ADJ_PETSC_DATA_STRUCTURES_H

#include <string.h>
#include "adj_constants.h"
#include "libadjoint/adj_data_structures.h"
#include "uthash.h"


/* void adj_set_petsc_data_callbacks(adj_adjointer adjointer); */

void petsc_vec_duplicate_proc(adj_vector x, adj_vector *newx);
void petsc_vec_axpy_proc(adj_vector *y, adj_scalar alpha, adj_vector x);
void petsc_vec_destroy_proc(adj_vector *x); 
void petsc_vec_setvalues_proc(adj_vector *vec, adj_scalar scalars[]); 
void petsc_vec_getsize_proc(adj_vector vec, int *sz);
void petsc_vec_divide_proc(adj_vector numerator, adj_vector denominator, adj_vector *output);

void petsc_mat_getvecs_proc(adj_matrix mat, adj_vector *left);
void petsc_mat_axpy_proc(adj_matrix *Y, adj_scalar alpha, adj_matrix X);
void petsc_mat_duplicate_proc(adj_matrix matin, adj_matrix *matout);
void petsc_mat_destroy_proc(adj_matrix *mat);

#endif
