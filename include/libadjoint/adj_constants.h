#ifndef ADJ_CONSTANTS_H
#define ADJ_CONSTANTS_H

#define ADJ_NAMELEN (size_t)255

/* values for the .type field of an adj_variable */
#define ADJ_FORWARD 1
#define ADJ_ADJOINT 2
#define ADJ_TLM 2

/* normal or auxiliary */
#define ADJ_NORMAL_VARIABLE 0
#define ADJ_AUXILIARY_VARIABLE 1

/* options for the adjointer */
#define ADJ_ACTIVITY 0
#define ADJ_ISP_ORDER 1

#define ADJ_ACTIVITY_ADJOINT 0
#define ADJ_ACTIVITY_NOTHING 1

#define ADJ_ISP_FIRSTORDER 0
#define ADJ_ISP_SECONDORDER 1

/* operator callbacks */
#define ADJ_NBLOCK_COLOURING_SZ_CB 0
#define ADJ_NBLOCK_COLOURING_CB 1
#define ADJ_NBLOCK_ACTION_CB 2
#define ADJ_NBLOCK_DERIVATIVE_ACTION_CB 3
#define ADJ_NBLOCK_DERIVATIVE_ASSEMBLY_CB 4
#define ADJ_BLOCK_ACTION 5
#define ADJ_BLOCK_ASSEMBLY 6

#define ADJ_VEC_DUPLICATE_CB 10
#define ADJ_VEC_AXPY_CB 11
#define ADJ_VEC_DESTROY_CB 12
#define ADJ_VEC_DIVIDE_CB 13
#define ADJ_VEC_SETVALUES_CB 14

#define ADJ_MAT_DUPLICATE_CB 20
#define ADJ_MAT_AXPY_CB 21
#define ADJ_MAT_DESTROY_CB 22
#define ADJ_MAT_GETVECS_CB 23
#endif
