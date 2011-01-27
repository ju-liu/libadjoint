#include <stdlib.h>
#include "adj_data_structures.h"
#include "adj_variable_lookup.h"
#include "adj_test_tools.h"
#include "adj_error_handling.h"

void test_adj_add_variable_data() 
{
  int ierr;

  adj_variable_hash *hash=NULL;
  adj_variable_data* data = (adj_variable_data*) malloc(sizeof(adj_variable_data));

  adj_variable a;
  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &a);

  data->equation = 19;

  ierr=adj_add_variable_data(hash, &a, data);
  adj_test_assert(ierr==ADJ_ERR_OK, "Should have worked");
  data = NULL;

  ierr=adj_find_variable_data(hash, &a, &data);
  adj_test_assert(ierr==ADJ_ERR_OK, "Should have worked");
  
/*  call adj_add_variable_data(adjointer, variable, data, ierr)
  call adj_test_assert(ierr == 0, "initial ierr")

  call adj_add_variable_data(adjointer, variable, data, ierr)
  call adj_test_assert(ierr /= 0, "second ierr")
  call adj_destroy_hash(adjointer) */
}
