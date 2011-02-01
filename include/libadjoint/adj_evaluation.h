#ifndef ADJ_EVALUATION_H
#define ADJ_EVALUATION_H

#include "libadjoint/adj_data_structures.h"
#include "libadjoint/adj_error_handling.h"
#include "libadjoint/adj_adjointer_routines.h"

int adj_evaluate_block_action(adj_adjointer* adjointer, adj_block block, adj_vector input, adj_vector* output);
int adj_evaluate_block_assembly(adj_adjointer* adjointer, adj_block block, adj_matrix *output, adj_vector* rhs);
int adj_evaluate_nonlinear_derivative_action(adj_adjointer* adjointer, adj_nonlinear_block_derivative* derivatives, adj_vector value, adj_vector rhs);

#endif
