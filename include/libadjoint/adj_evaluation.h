#ifndef ADJ_EVALUATION_H
#define ADJ_EVALUATION_H

#include "adj_data_structures.h"
#include "adj_error_handling.h"
#include "adj_adjointer_routines.h"
#include "adj_debug.h"

#ifndef ADJ_HIDE_FROM_USER
int adj_evaluate_block_action(adj_adjointer* adjointer, adj_block block, adj_vector input, adj_vector* output);
int adj_evaluate_block_assembly(adj_adjointer* adjointer, adj_block block, adj_matrix *output, adj_vector* rhs);
int adj_evaluate_nonlinear_action(adj_adjointer* adjointer, void (*nonlinear_action_func)(int nvar, adj_variable* variables, adj_vector* dependencies,
     adj_vector input, void* context, adj_vector* output),adj_nonlinear_block nonlinear_block, adj_vector input, adj_variable* perturbed_var,
     adj_vector* perturbation, adj_vector* output);
int adj_evaluate_nonlinear_derivative_action(adj_adjointer* adjointer, int nderivatives, adj_nonlinear_block_derivative* derivatives, adj_vector value, adj_vector model_rhs, adj_vector* rhs);
int adj_evaluate_nonlinear_derivative_action_supplied(adj_adjointer* adjointer, void (*nonlinear_derivative_action_func)(int nvar, adj_variable* variables, 
     adj_vector* dependencies, adj_variable derivative, adj_vector contraction, int hermitian, adj_vector input, adj_scalar coefficient, void* context, adj_vector* output),
     adj_nonlinear_block_derivative derivative, adj_vector value, adj_vector* rhs);
int adj_evaluate_nonlinear_derivative_action_isp(adj_adjointer* adjointer, void (*nonlinear_action_func)(int nvar, adj_variable*
     variables, adj_vector* dependencies, adj_vector input, void* context, adj_vector* output), adj_nonlinear_block_derivative derivative, 
     adj_vector value, adj_vector* rhs);
int adj_evaluate_nonlinear_colouring(adj_adjointer* adjointer, adj_nonlinear_block nonlinear_block, adj_variable derivative,
    void (*nonlinear_colouring_func)(int nvar, adj_variable* variables, adj_vector* dependencies, adj_variable derivative, void* context, int sz, int* colouring),
    int sz, int* colouring);
int adj_evaluate_functional(adj_adjointer* adjointer, int timestep, char* functional, adj_scalar* output);
int adj_evaluate_functional_derivative(adj_adjointer* adjointer, adj_variable variable, char* functional, adj_vector* output, int* has_output);
int adj_evaluate_forward_source(adj_adjointer* adjointer, int equation, adj_vector* output, int* has_output);
#endif

#endif
