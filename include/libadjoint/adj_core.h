#include "adj_data_structures.h"
#include "adj_error_handling.h"
#include "adj_variable_lookup.h"
#include "adj_adjointer_routines.h"
#include "adj_evaluation.h"
#include "adj_simplification.h"

int adj_get_adjoint_equation(adj_adjointer* adjointer, int equation, int functional, adj_matrix* lhs, adj_vector* rhs, adj_variable* var);
