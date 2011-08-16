#include "adj_data_structures.h"
#include "adj_error_handling.h"
#include "adj_variable_lookup.h"
#include "adj_adjointer_routines.h"
#include "adj_evaluation.h"
#include "adj_simplification.h"
#include "revolve_c.h"

int adj_get_adjoint_equation(adj_adjointer* adjointer, int equation, char* functional, adj_matrix* lhs, adj_vector* rhs, adj_variable* var);
int adj_get_adjoint_solution(adj_adjointer* adjointer, int equation, char* functional, adj_vector* soln, adj_variable* var);
int adj_get_forward_equation(adj_adjointer* adjointer, int equation, adj_matrix* lhs, adj_vector* rhs, adj_variable* var);
