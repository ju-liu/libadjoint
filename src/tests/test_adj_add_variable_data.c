#include <stdlib.h>
#include "libadjoint/adj_data_structures.h"
#include "libadjoint/adj_variable_lookup.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_error_handling.h"
#include "libadjoint/adj_test_main.h"

void test_adj_add_variable_data(void)
{
  int ierr;

  adj_variable_hash *hash=NULL;
  adj_variable_data* data = (adj_variable_data*) malloc(sizeof(adj_variable_data));
  data->equation = 19;

  adj_variable a;
  adj_create_variable("Velocity", 0, 0, ADJ_NORMAL_VARIABLE, &a);

  ierr=adj_add_variable_data(&hash, &a, data);
  adj_test_assert(ierr==ADJ_OK, "Should have worked");
  data = NULL;

  ierr=adj_find_variable_data(&hash, &a, &data);
  adj_test_assert(ierr==ADJ_OK, "Should have worked");
  adj_test_assert(data->equation == 19, "Should have worked");

  ierr=adj_add_variable_data(&hash, &a, data);
  adj_test_assert(ierr!=ADJ_OK, "Should not have worked");
  
  adj_destroy_hash(&hash);
  free(data);
}
