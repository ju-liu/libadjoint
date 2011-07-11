#include "adj_data_structures.h"
#include "adj_adjointer_routines.h"
#include "adj_error_handling.h"

#ifndef ADJ_HIDE_FROM_USER
int adj_simplify_derivatives(adj_adjointer* adjointer, int ninput, adj_nonlinear_block_derivative* input, int* noutput, adj_nonlinear_block_derivative** output);
int adj_simplification_compare(adj_nonlinear_block_derivative d1, adj_nonlinear_block_derivative d2);
void adj_simplification_merge(adj_adjointer* adjointer, adj_nonlinear_block_derivative* d1, adj_nonlinear_block_derivative* d2, int merged);
#endif
