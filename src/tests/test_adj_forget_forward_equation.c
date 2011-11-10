#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_test_main.h"
#include "libadjoint/adj_adjointer_visualisation.h"
#include "libadjoint/adj_core.h"

#ifndef HAVE_PETSC
void test_adj_forget_forward_equation(void)
{
  adj_test_assert(1 == 1, "Don't have PETSc so can't run this test.");
}
#else
#include "libadjoint/adj_petsc_data_structures.h"

void burgers_operator_assembly(int ndepends, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_scalar coefficient, void* context, adj_matrix* output, adj_vector* rhs);
void timestepping_operator_action(int ndepends, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_scalar coefficient, adj_vector input, void* context, adj_vector* output);
void advection_derivative_operator_action(int ndepends, adj_variable* variables, adj_vector* dependencies, adj_variable derivative, adj_vector contraction, int hermitian, adj_vector input, adj_scalar coefficient, void* context, adj_vector* output);
void source(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, void* context, adj_vector* output, int* has_output);

void test_adj_forget_forward_equation(void)
{
  int steps = 5;
  int snaps = 4;
  int snaps_in_ram = 1;
  int timestep, nb_eqs;
  int cs;
  adj_adjointer adjointer;
  adj_nonlinear_block V;
  adj_block B[2], I;
  adj_variable u[2];
  adj_equation eqn;
  int ierr;
  Vec vec;
  adj_vector value;
  int dim=2;
  adj_storage_data storage;
  char filename_fwd[ADJ_NAME_LEN];
  int checkpoint_strategy[steps];
  adj_variable_data *var_data;

  ierr = adj_create_adjointer(&adjointer);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  ierr = adj_set_checkpoint_strategy(&adjointer, ADJ_CHECKPOINT_REVOLVE_OFFLINE);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  ierr = adj_set_revolve_options(&adjointer, steps, snaps, snaps_in_ram, ADJ_TRUE);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  /* Register callbacks */
  ierr = adj_set_petsc_data_callbacks(&adjointer);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_register_forward_source_callback(&adjointer, source); 
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr =  adj_register_operator_callback(&adjointer, ADJ_BLOCK_ASSEMBLY_CB, "BurgersOperator", (void (*)(void)) burgers_operator_assembly);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr =  adj_register_operator_callback(&adjointer, ADJ_BLOCK_ASSEMBLY_CB, "IdentityOperator", (void (*)(void)) burgers_operator_assembly);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr =  adj_register_operator_callback(&adjointer, ADJ_BLOCK_ACTION_CB, "TimesteppingOperator", (void (*)(void)) timestepping_operator_action);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr =  adj_register_operator_callback(&adjointer, ADJ_NBLOCK_DERIVATIVE_ACTION_CB, "AdvectionOperator", (void (*)(void)) advection_derivative_operator_action);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");


  /* Our solution */
  VecCreateSeq(PETSC_COMM_SELF, dim, &vec);
  VecSet(vec, 1.0);
  value = petsc_vec_to_adj_vector(&vec);

  /* Initial condition */
  timestep = 0;
  adj_create_variable("Velocity", timestep, 0, ADJ_NORMAL_VARIABLE, &u[1]);
  adj_create_block("IdentityOperator", NULL, NULL, &I);
  ierr = adj_create_equation(u[1], 1, &I, &u[1], &eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_register_equation(&adjointer, eqn, &cs);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  checkpoint_strategy[0] = cs;
  adj_test_assert(cs!=ADJ_CHECKPOINT_STORAGE_MEMORY, "Offline checkpointing should not use storage_memory");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  /* Annotate the model */
  for (timestep=1;timestep<steps; timestep++)
  {
    u[0] = u[1];
    adj_create_variable("Velocity", timestep, 0, ADJ_NORMAL_VARIABLE, &u[1]);

    adj_create_nonlinear_block("AdvectionOperator", 1, &u[0], NULL, &V);
    adj_nonlinear_block_set_coefficient(&V, 0.5);
    adj_create_block("TimesteppingOperator", &V, NULL, &B[0]);
    adj_create_block("BurgersOperator", &V, NULL, &B[1]);
    ierr = adj_create_equation(u[1], 2, B, u, &eqn);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");
    ierr = adj_register_equation(&adjointer, eqn, &cs);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");
    checkpoint_strategy[timestep] = cs;
    ierr = adj_destroy_equation(&eqn);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");
    adj_destroy_nonlinear_block(&V);
    adj_destroy_block(&B[0]);
    adj_destroy_block(&B[1]);
  }


  ierr = adj_equation_count(&adjointer, &nb_eqs);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  adj_test_assert(nb_eqs  == steps, "Number of timesteps should be the same the number of registered equations");


  /* Replay this annotation */
  for (timestep=0;timestep<steps; timestep++)
  {
    printf("Checkpoint strategy for timestep %i: %i\n", timestep, checkpoint_strategy[timestep]);
    ierr = adj_get_forward_solution(&adjointer, timestep, &value, &u[0]);
    adj_chkierr(ierr);

    if (checkpoint_strategy[timestep]==ADJ_CHECKPOINT_STORAGE_DISK)
    {
      ierr = adj_storage_disk(value, &storage);
      adj_test_assert(ierr == ADJ_OK, "Should have worked");
      /* A checkpoint needs the the old velocity only */
      ierr = adj_record_variable(&adjointer, u[0], storage);
      adj_test_assert(ierr == ADJ_OK, "Should have worked");
    }
     ierr = adj_storage_memory_copy(value, &storage);
     adj_test_assert(ierr == ADJ_OK, "Should have worked");
    ierr = adj_record_variable(&adjointer, u[0], storage);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");

    /* Output */
    sprintf(filename_fwd, "test_adj_forget_forward_equation_%i_before_forget.html", timestep);
    ierr = adj_adjointer_to_html(&adjointer, filename_fwd, ADJ_FORWARD);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");

    adj_forget_forward_equation(&adjointer, timestep);

    /* Output after forget */
    sprintf(filename_fwd, "test_adj_forget_forward_equation_%i_after_forget.html", timestep);
    ierr = adj_adjointer_to_html(&adjointer, filename_fwd, ADJ_FORWARD);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");
  }

  /* Besides checkpoints, no variables should be stored anymore ... */
  var_data=adjointer.vardata.firstnode;
  while(var_data!=NULL)
  {
     adj_test_assert(var_data->storage.storage_memory_has_value==var_data->storage.storage_memory_is_checkpoint, "Should have worked");
     adj_test_assert(var_data->storage.storage_disk_has_value==var_data->storage.storage_disk_is_checkpoint, "Should have worked");

     var_data=var_data->next;
  }

}

void source(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, void* context, adj_vector* output, int* has_output)
{
  (void) adjointer;
  (void) variable;
  (void) ndepends;
  (void) variables;
  (void) dependencies;
  (void) context;
  (void) output;
  *has_output=ADJ_FALSE;
}

void burgers_operator_assembly(int ndepends, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_scalar coefficient, void* context, adj_matrix* output_mat, adj_vector* rhs_vec)
{
  int dim=2;
  Mat *output;
  Vec *rhs, *ones;
  (void) ndepends;
  (void) variables;
  (void) dependencies;
  (void) hermitian;
  (void) coefficient;
  (void) context;

  output = (Mat*) malloc(sizeof(Mat));
  rhs = (Vec*) malloc(sizeof(Vec));
  ones = (Vec*) malloc(sizeof(Vec));

  MatCreateSeqAIJ(PETSC_COMM_SELF, dim, dim, 1, PETSC_NULL, output);
  VecCreateSeq(PETSC_COMM_SELF, dim, ones);
  VecCreateSeq(PETSC_COMM_SELF, dim, rhs);
  VecZeroEntries(*rhs);
  VecSet(*ones, 1.0);
  MatDiagonalSet(*output, *ones, INSERT_VALUES);
  VecDestroy(*ones);
  if (hermitian == ADJ_TRUE) 
    MatHermitianTranspose(*output, MAT_REUSE_MATRIX, output);
  *rhs_vec = petsc_vec_to_adj_vector(rhs);
  *output_mat = petsc_mat_to_adj_matrix(output);
}

void timestepping_operator_action(int ndepends, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_scalar coefficient, adj_vector input, void* context, adj_vector* output)
{
  Vec *output_vec;
  (void) ndepends;
  (void) variables;
  (void) dependencies;
  (void) hermitian;
  (void) coefficient;
  (void) context;

  output_vec = (Vec*) malloc(sizeof(Vec));

  VecDuplicate(petsc_vec_from_adj_vector(input), output_vec);
  VecCopy(petsc_vec_from_adj_vector(input), *output_vec);
  *output = petsc_vec_to_adj_vector(output_vec);
}

void advection_derivative_operator_action(int ndepends, adj_variable* variables, adj_vector* dependencies, adj_variable derivative, adj_vector contraction, int hermitian, adj_vector input, adj_scalar coefficient, void* context, adj_vector* output)
{
  Vec *output_vec;
  (void) ndepends;
  (void) variables;
  (void) dependencies;
  (void) derivative; 
  (void) contraction;
  (void) hermitian;
  (void) coefficient;
  (void) context;

  output_vec = (Vec*) malloc(sizeof(Vec));

  VecDuplicate(petsc_vec_from_adj_vector(input), output_vec);
  VecCopy(petsc_vec_from_adj_vector(input), *output_vec);
  *output = petsc_vec_to_adj_vector(output_vec);
}
#endif
