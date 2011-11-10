#ifndef ADJ_PETSC_DATA_STRUCTURES_H
#define ADJ_PETSC_DATA_STRUCTURES_H

#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include "adj_constants.h"
#include "adj_adjointer_routines.h"
#include "adj_data_structures.h"
#include "uthash.h"
#ifdef HAVE_PETSC
#include "libadjoint/adj_petsc.h"
#endif

int adj_set_petsc_data_callbacks(adj_adjointer *adjointer);

void petsc_vec_duplicate_proc(adj_vector x, adj_vector *newx);
void petsc_vec_axpy_proc(adj_vector *y, adj_scalar alpha, adj_vector x);
void petsc_vec_destroy_proc(adj_vector *x); 
void petsc_vec_setvalues_proc(adj_vector *vec, adj_scalar scalars[]); 
void petsc_vec_getsize_proc(adj_vector vec, int *sz);
void petsc_vec_divide_proc(adj_vector *numerator, adj_vector denominator);
void petsc_vec_getnorm_proc(adj_vector vec, adj_scalar* norm);
void petsc_vec_set_random_proc(adj_vector* x);
void petsc_vec_dot_product_proc(adj_vector x, adj_vector y, adj_scalar* val);
void petsc_vec_to_file_proc(adj_vector x, char* filename); 
void petsc_vec_from_file_proc(adj_vector* x, char* filename);

void petsc_mat_getvec_proc(adj_matrix mat, adj_vector *left);
void petsc_mat_axpy_proc(adj_matrix *Y, adj_scalar alpha, adj_matrix X);
void petsc_mat_duplicate_proc(adj_matrix matin, adj_matrix *matout);
void petsc_mat_destroy_proc(adj_matrix *mat);
void petsc_solve_proc(adj_variable var, adj_matrix mat, adj_vector rhs, adj_vector *soln); 

#ifdef HAVE_PETSC
adj_vector petsc_vec_to_adj_vector(Vec* v);
Vec petsc_vec_from_adj_vector(adj_vector vector);
adj_matrix petsc_mat_to_adj_matrix(Mat* v);
Mat petsc_mat_from_adj_matrix(adj_matrix matrix);
#endif

#endif
