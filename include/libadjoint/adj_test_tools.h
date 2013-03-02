#ifndef ADJ_TEST_TOOLS_H
#define ADJ_TEST_TOOLS_H
#include "adj_data_structures.h"

#ifdef __cplusplus
extern "C" {
#endif

void adj_test_assert(int passed, char *testdesc);
int adj_sizeof_adjointer(void);

#ifdef __cplusplus
}
#endif
#endif

