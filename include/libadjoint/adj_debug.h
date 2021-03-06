#ifndef ADJ_DEBUG_H
#define ADJ_DEBUG_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#undef I
#include "adj_data_structures.h"
#include "adj_error_handling.h"
#include "adj_evaluation.h"

#ifdef __cplusplus
#include <complex>
#define COMPLEX_ABS(x) (abs( (complex<double>) (x)))
#else
#include <complex.h>
#define COMPLEX_ABS(x) (cabs( (double complex) (x)))
#endif

int adj_adjointer_check_consistency(adj_adjointer* adjointer);
int adj_adjointer_check_checkpoints(adj_adjointer* adjointer); /* Checks that all required variables for the adjoint solve with revolve are stored as checkpoints. Should be called just before the adjoint main loop */
#ifndef ADJ_HIDE_FROM_USER
int adj_test_block_action_transpose(adj_adjointer* adjointer, adj_block block, adj_vector model_input, adj_vector model_output, int N, adj_scalar tol);
int adj_test_nonlinear_derivative_action_transpose(adj_adjointer* adjointer, adj_nonlinear_block_derivative nonlinear_block_derivative, adj_vector model_input, adj_vector model_output, int N, adj_scalar tol);
int adj_test_nonlinear_derivative_action_consistency(adj_adjointer* adjointer, adj_nonlinear_block_derivative nonlinear_block_derivative, adj_variable deriv_var, int N);
#endif

#endif
