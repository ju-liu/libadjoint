#ifndef ADJ_CONSTANTS_H
#define ADJ_CONSTANTS_H

/* You cannot use comments of more than one line in this file */
/* Maximum length of a name for something in libadjoint */
#define ADJ_NAME_LEN 255
/* Used in the string-to-string hash table */
#define ADJ_DICT_LEN 8192

#define adj_scalar double
#define adj_scalar_f real(kind=c_double)
#define adj_scalar_f_cast(X) real((X), kind=c_double)
#define ADJ_SCALAR_EPS 1.0e-13

#define ADJ_TRUE 1
#define ADJ_FALSE 0

/* values for the .type field of an adj_variable */
#define ADJ_FORWARD 1
#define ADJ_ADJOINT 2 /* adjoint variable */
#define ADJ_TLM 3 /* tangent linear variable */
#define ADJ_SOA 4 /* second order adjoint */

/* normal or auxiliary */
#define ADJ_NORMAL_VARIABLE 0
#define ADJ_AUXILIARY_VARIABLE 1

/* options for the adjointer */
#define ADJ_NO_OPTIONS 3
#define ADJ_ACTIVITY 0
#define ADJ_ISP_ORDER 1
#define ADJ_CHECKPOINT_STRATEGY 2

/* whichever value is zero defines the default */
#define ADJ_ACTIVITY_ADJOINT 0
#define ADJ_ACTIVITY_NOTHING 1

#define ADJ_CHECKPOINT_NONE 0
#define ADJ_CHECKPOINT_REVOLVE_OFFLINE 1
#define ADJ_CHECKPOINT_REVOLVE_MULTISTAGE 2
#define ADJ_CHECKPOINT_REVOLVE_ONLINE 3

#define ADJ_ISP_SECONDORDER 0
#define ADJ_ISP_FIRSTORDER 1

/* storage strategies for checkpointing */
#define ADJ_CHECKPOINT_STORAGE_NONE 0
#define ADJ_CHECKPOINT_STORAGE_MEMORY 1
#define ADJ_CHECKPOINT_STORAGE_DISK 2

/* storage strategies */
#define ADJ_STORAGE_MEMORY_COPY 0
#define ADJ_STORAGE_MEMORY_INCREF 1

/* operator callbacks */
#define ADJ_NBLOCK_COLOURING_CB 1
#define ADJ_NBLOCK_ACTION_CB 2
#define ADJ_NBLOCK_DERIVATIVE_ACTION_CB 3
#define ADJ_NBLOCK_DERIVATIVE_ASSEMBLY_CB 4
#define ADJ_BLOCK_ACTION_CB 5
#define ADJ_BLOCK_ASSEMBLY_CB 6
/* if you add a new one, you must update the table in src/adj_adjointer_routines.c */

#define ADJ_VEC_DUPLICATE_CB 10
#define ADJ_VEC_AXPY_CB 11
#define ADJ_VEC_DESTROY_CB 12
#define ADJ_VEC_DIVIDE_CB 13
#define ADJ_VEC_SET_VALUES_CB 14
#define ADJ_VEC_GET_VALUES_CB 15
#define ADJ_VEC_GET_SIZE_CB 16
#define ADJ_VEC_GET_NORM_CB 17
#define ADJ_VEC_DOT_PRODUCT_CB 18
#define ADJ_VEC_SET_RANDOM_CB 19
#define ADJ_VEC_WRITE_CB 20
#define ADJ_VEC_READ_CB 21
#define ADJ_VEC_DELETE_CB 22

#define ADJ_MAT_DUPLICATE_CB 30
#define ADJ_MAT_AXPY_CB 31
#define ADJ_MAT_DESTROY_CB 32
#define ADJ_MAT_ACTION_CB 33

#define ADJ_SOLVE_CB 40

/* prealloc constant */
#define ADJ_PREALLOC_SIZE 1

/* value for unset variables */
#define ADJ_UNSET -666
#endif
