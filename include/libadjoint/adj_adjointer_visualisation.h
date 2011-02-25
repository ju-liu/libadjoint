#ifndef ADJ_ADJOINTER_VISUALISATION_H
#define ADJ_ADJOINTER_VISUALISATION_H

#include <stdio.h>
#include <stdlib.h>
#include "adj_data_structures.h"
#include "adj_constants.h"
#include "adj_error_handling.h"
#include "adj_variable_lookup.h"

int adj_adjointer_to_html(adj_adjointer* adjointer, char* filename, int type);

#endif
