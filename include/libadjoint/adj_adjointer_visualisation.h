#ifndef ADJ_ADJOINTER_VISUALISATION_H
#define ADJ_ADJOINTER_VISUALISATION_H

#include <stdio.h>
#include <stdlib.h>
#include "adj_data_structures.h"
#include "adj_constants.h"
#include "adj_error_handling.h"
#include "adj_variable_lookup.h"

int adj_adjointer_to_html(adj_adjointer* adjointer, char* filename, int type);

#ifndef ADJ_HIDE_FROM_USER
void adj_html_css(FILE* fp);
void adj_html_header(FILE *fp);
void adj_html_footer(FILE *fp);
void adj_html_table_begin(FILE* fp);
void adj_html_table_end(FILE* fp);
void adj_html_write_row(FILE* fp, char** strings, char** desc, int nb_strings);
int adj_html_find_column_index(adj_adjointer* adjointer, adj_variable* variable, int* col);
void adj_html_vars(FILE* fp, adj_adjointer* adjointer, int type);
int adj_html_eqn(FILE* fp, adj_adjointer* adjointer, adj_equation adj_eqn);
int adj_html_adjoint_eqn(FILE* fp, adj_adjointer* adjointer, adj_equation fwd_eqn);
int adj_html_adjoint_system(adj_adjointer* adjointer, char* filename);
int adj_html_forward_system(adj_adjointer* adjointer, char* filename);
#endif

#endif