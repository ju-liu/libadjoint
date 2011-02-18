#ifndef ADJ_EVALUATION_H
#define ADJ_EVALUATION_H

#include "adj_data_structures.h"
#include "adj_error_handling.h"
#include "adj_adjointer_routines.h"


#ifndef ADJ_HIDE_FROM_USER
int adj_evaluate_block_action(adj_adjointer* adjointer, adj_block block, adj_vector input, adj_vector* output);
int adj_evaluate_block_assembly(adj_adjointer* adjointer, adj_block block, adj_matrix *output, adj_vector* rhs);
int adj_evaluate_nonlinear_derivative_action(adj_adjointer* adjointer, adj_nonlinear_block_derivative* derivatives, adj_vector value, adj_vector rhs);
int adj_evaluate_functional(adj_adjointer* adjointer, adj_variable variable, char* functional, adj_vector* output);
#endif

#endif
