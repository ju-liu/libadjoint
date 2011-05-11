#ifndef ADJ_DEBUG_H
#define ADJ_DEBUG_H

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#undef I
#include "adj_data_structures.h"
#include "adj_error_handling.h"
#include "adj_evaluation.h"

int adj_adjointer_check_consistency(adj_adjointer* adjointer);
int adj_test_action_transpose(adj_adjointer* adjointer, adj_block block, adj_vector model_input, adj_vector model_output, int N, adj_scalar tol);

#endif
