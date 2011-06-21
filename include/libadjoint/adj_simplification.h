#include "adj_data_structures.h"
#include "adj_adjointer_routines.h"

#ifndef ADJ_HIDE_FROM_USER
void adj_simplify_derivatives(adj_adjointer* adjointer, adj_nonlinear_block_derivative* input, adj_nonlinear_block_derivative** output);
int adj_simplification_compare(adj_nonlinear_block_derivative d1, adj_nonlinear_block_derivative d2);
void adj_simplficiation_merge(adj_adjointer* adjointer, adj_nonlinear_block_derivative* d1, adj_nonlinear_block_derivative d2, int merged);
#endif
