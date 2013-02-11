#include "adj_data_structures.h"
#include "adj_error_handling.h"
#include "adj_variable_lookup.h"
#include "adj_adjointer_routines.h"
#include "adj_evaluation.h"
#include "adj_simplification.h"
#include "revolve_c.h"

int adj_get_adjoint_equation(adj_adjointer* adjointer, int equation, char* functional, adj_matrix* lhs, adj_vector* rhs, adj_variable* adj_var);
int adj_get_adjoint_solution(adj_adjointer* adjointer, int equation, char* functional, adj_vector* soln, adj_variable* adj_var);
int adj_get_forward_equation(adj_adjointer* adjointer, int equation, adj_matrix* lhs, adj_vector* rhs, adj_variable* fwd_var);
int adj_get_forward_solution(adj_adjointer* adjointer, int equation, adj_vector* soln, adj_variable* fwd_var);
int adj_get_tlm_equation    (adj_adjointer* adjointer, int equation, char* parameter,  adj_matrix* lhs, adj_vector* rhs, adj_variable* tlm_var);
int adj_get_tlm_solution    (adj_adjointer* adjointer, int equation, char* parameter,  adj_vector* soln, adj_variable* tlm_var);
int adj_get_soa_equation    (adj_adjointer* adjointer, int equation, char* functional, char* parameter,  adj_matrix* lhs, adj_vector* rhs, adj_variable* soa_var);
int adj_get_soa_solution    (adj_adjointer* adjointer, int equation, char* functional, char* parameter,  adj_vector* soln, adj_variable* soa_var);

#ifndef ADJ_HIDE_FROM_USER
int adj_replay_forward_equations(adj_adjointer* adjointer, int start_equation, int stop_equation, int checkpoint_last_timestep);
int adj_revolve_to_adjoint_equation(adj_adjointer* adjointer, int equation);
#endif
