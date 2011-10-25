#include "libadjoint/adj_adjointer_routines.h"
#include "libadjoint/adj_test_tools.h"
#include "libadjoint/adj_test_main.h"
#include "libadjoint/adj_adjointer_visualisation.h"
#include "libadjoint/adj_core.h"

#ifndef HAVE_PETSC
void test_checkpoint_revolve_online(void)
{
  adj_test_assert(1 == 1, "Don't have PETSc so can't run this test.");
}
#else
#include "libadjoint/adj_petsc_data_structures.h"

void drag_derivative(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output);
void burgers_operator_assembly(int ndepends, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_scalar coefficient, void* context, adj_matrix* output, adj_vector* rhs);
void timestepping_operator_action(int ndepends, adj_variable* variables, adj_vector* dependencies, int hermitian, adj_scalar coefficient, adj_vector input, void* context, adj_vector* output);
void advection_derivative_operator_action(int ndepends, adj_variable* variables, adj_vector* dependencies, adj_variable derivative, adj_vector contraction, int hermitian, adj_vector input, adj_scalar coefficient, void* context, adj_vector* output);
void source(adj_adjointer* adjointer, adj_variable variable, int ndepends, adj_variable* variables, adj_vector* dependencies, void* context, adj_vector* output, int* has_output);

int test_checkpoints(adj_adjointer *adjointer, int nb_expected_vars, char expected_vars[][ADJ_NAME_LEN], int* memory_has_value, int* memory_is_checkpoint, int* disk_has_value, int* disk_is_checkpoint);
void get_expected_values(int timestep, int* nb_expected_vars, char expected_vars[][ADJ_NAME_LEN], int* memory_has_value, int* memory_is_checkpoint, int* disk_has_value, int* disk_is_checkpoint);

void test_checkpoint_revolve_online(void)
{
  int steps = 20;
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
  char filename_adj[ADJ_NAME_LEN];

  /* Expected results */
  int nb_expected_vars;
  char expected_vars[steps][ADJ_NAME_LEN];
  int memory_has_value[steps];
  int memory_is_checkpoint[steps];
  int disk_has_value[steps];
  int disk_is_checkpoint[steps];

  /* Variables for the adjoint */  
  adj_variable lambda;

  ierr = adj_create_adjointer(&adjointer);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  ierr = adj_set_checkpoint_strategy(&adjointer, ADJ_CHECKPOINT_REVOLVE_ONLINE);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  ierr = adj_set_revolve_options(&adjointer, -1, snaps, snaps_in_ram, ADJ_TRUE); /* We do not tell revolve how many steps we are going to do */
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  /* Register callbacks */
  ierr = adj_set_petsc_data_callbacks(&adjointer);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_register_forward_source_callback(&adjointer, source); 
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_register_functional_derivative_callback(&adjointer, "Drag", drag_derivative); 
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
  printf("Checkpoint strategy for timestep %i: %i\n", 0, cs);
  adj_test_assert(cs!=ADJ_CHECKPOINT_STORAGE_MEMORY, "Offline checkpointing should not use storage_memory");
  ierr = adj_destroy_equation(&eqn);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  if (cs==ADJ_CHECKPOINT_STORAGE_DISK)
  {
  	/* We do not need anything for the initial checkpoint */
  }

  ierr = adj_storage_memory_copy(value, &storage);
	adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_record_variable(&adjointer, u[1], storage);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  /* A typical time loop */
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
    printf("Checkpoint strategy for timestep %i: %i\n", timestep, cs);
    ierr = adj_destroy_equation(&eqn);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");
    adj_destroy_nonlinear_block(&V);
    adj_destroy_block(&B[0]);
    adj_destroy_block(&B[1]);

    /* We need a checkpoint at the last timestep */
    if (timestep==steps-1)
    	cs=ADJ_CHECKPOINT_STORAGE_MEMORY;

    if (cs==ADJ_CHECKPOINT_STORAGE_DISK)
    {
  		ierr = adj_storage_disk(value, &storage);
  		adj_test_assert(ierr == ADJ_OK, "Should have worked");
  		ierr = adj_storage_set_checkpoint(&storage, ADJ_TRUE);
  		adj_test_assert(ierr == ADJ_OK, "Should have worked");
  		/* A checkpoint needs the the old velocity only */
  		ierr = adj_record_variable(&adjointer, u[0], storage);
  		adj_test_assert(ierr == ADJ_OK, "Should have worked");
    }
    else if (cs==ADJ_CHECKPOINT_STORAGE_MEMORY)
    {
			ierr = adj_storage_memory_copy(value, &storage);
			adj_test_assert(ierr == ADJ_OK, "Should have worked");
			ierr = adj_storage_set_checkpoint(&storage, ADJ_TRUE);
			adj_test_assert(ierr == ADJ_OK, "Should have worked");
			ierr = adj_record_variable(&adjointer, u[0], storage);
			adj_test_assert(ierr == ADJ_OK, "Should have worked");
    }

    /* solve for u[1] .... */

    /* We need to tell revolve how to save results in memory while replaying */
		ierr = adj_set_storage_memory_copy(&adjointer, &u[1]);
		adj_test_assert(ierr == ADJ_OK, "Should have worked");
  }

  /* We can record the solution of the last timestep, but we do not have to: revolve will do it for us */
  /*{
		ierr = adj_storage_memory_copy(value, &storage);
		adj_test_assert(ierr == ADJ_OK, "Should have worked");
		ierr = adj_storage_set_checkpoint(&storage, ADJ_TRUE);
		adj_test_assert(ierr == ADJ_OK, "Should have worked");
		ierr = adj_record_variable(&adjointer, u[1], storage);
		adj_test_assert(ierr == ADJ_OK, "Should have worked");
  }*/

  ierr = adj_equation_count(&adjointer, &nb_eqs);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  adj_test_assert(nb_eqs  == steps, "Number of timesteps should be the same the number of registered equations");

  /* Forget the variables of the previous last timestep */
  ierr = adj_forget_forward_equation(&adjointer, steps-2);
  //ierr = adj_forget_forward_equation_until(&adjointer, 6, 7);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");


  ierr = adj_adjointer_to_html(&adjointer, "test_revolve_checkpoint_forward.html", ADJ_FORWARD);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");
  ierr = adj_adjointer_to_html(&adjointer, "test_revolve_checkpoint_adjoint.html", ADJ_ADJOINT);
  adj_test_assert(ierr == ADJ_OK, "Should have worked");

  ierr = adj_adjointer_check_checkpoints(&adjointer);
  adj_chkierr(ierr);

  /* A typical adjoint solve */
  for (timestep=steps-1;timestep>=0; timestep--)
  {
    ierr = adj_get_adjoint_solution(&adjointer, timestep, "Drag", &value, &lambda);
    adj_chkierr(ierr);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");

    ierr = adj_storage_memory_copy(value, &storage);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");

    ierr = adj_record_variable(&adjointer, lambda, storage);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");

    /* Now forget unnecessary variables */
    ierr = adj_forget_adjoint_equation(&adjointer, timestep);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");

    sprintf(filename_fwd, "test_revolve_checkpoint_forward_%i.html", timestep);
    sprintf(filename_adj, "test_revolve_checkpoint_adjoint_%i.html", timestep);

    ierr = adj_adjointer_to_html(&adjointer, filename_fwd, ADJ_FORWARD);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");
    ierr = adj_adjointer_to_html(&adjointer, filename_adj, ADJ_ADJOINT);
    adj_test_assert(ierr == ADJ_OK, "Should have worked");

    /* Compare the recorded forward variables with the expected results */
    if ((timestep>=14) || (timestep<=5))
    {
			/*get_expected_values(timestep, &nb_expected_vars, expected_vars, memory_has_value, memory_is_checkpoint, disk_has_value, disk_is_checkpoint);
			ierr = test_checkpoints(&adjointer, nb_expected_vars, expected_vars, memory_has_value, memory_is_checkpoint, disk_has_value, disk_is_checkpoint);
			if (ierr!=ADJ_OK)
				printf("Error in timestep %i", timestep);
			adj_test_assert(ierr == ADJ_OK, "Should have worked");
			*/
    }
  }

}

void get_expected_values(int timestep, int* nb_expected_vars, char expected_vars[][ADJ_NAME_LEN], int* memory_has_value, int* memory_is_checkpoint, int* disk_has_value, int* disk_is_checkpoint)
{
	switch(timestep)
	{
	case 19:
		*nb_expected_vars=4;
		strncpy(expected_vars[0], "Velocity:9:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[1], "Velocity:15:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[2], "Velocity:18:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[3], "Velocity:19:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_has_value[1]=0;
		memory_has_value[2]=1;
		memory_has_value[3]=1;
		memory_is_checkpoint[0]=0;
		memory_is_checkpoint[1]=0;
		memory_is_checkpoint[2]=1;
		memory_is_checkpoint[3]=1;
		disk_has_value[0]=1;
		disk_has_value[1]=1;
		disk_has_value[2]=0;
		disk_has_value[3]=0;
		disk_is_checkpoint[0]=1;
		disk_is_checkpoint[1]=1;
		disk_is_checkpoint[2]=0;
		disk_is_checkpoint[3]=0;
		break;
	case 18:
		*nb_expected_vars=4;
		strncpy(expected_vars[0], "Velocity:9:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[1], "Velocity:15:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[2], "Velocity:17:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[3], "Velocity:18:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_has_value[1]=0;
		memory_has_value[2]=1;
		memory_has_value[3]=1;
		memory_is_checkpoint[0]=0;
		memory_is_checkpoint[1]=0;
		memory_is_checkpoint[2]=1;
		memory_is_checkpoint[3]=1;
		disk_has_value[0]=1;
		disk_has_value[1]=1;
		disk_has_value[2]=0;
		disk_has_value[3]=0;
		disk_is_checkpoint[0]=1;
		disk_is_checkpoint[1]=1;
		disk_is_checkpoint[2]=0;
		disk_is_checkpoint[3]=0;
		break;
	case 17:
		*nb_expected_vars=4;
		strncpy(expected_vars[0], "Velocity:9:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[1], "Velocity:15:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[2], "Velocity:16:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[3], "Velocity:17:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_has_value[1]=0;
		memory_has_value[2]=1;
		memory_has_value[3]=1;
		memory_is_checkpoint[0]=0;
		memory_is_checkpoint[1]=0;
		memory_is_checkpoint[2]=1;
		memory_is_checkpoint[3]=1;
		disk_has_value[0]=1;
		disk_has_value[1]=1;
		disk_has_value[2]=0;
		disk_has_value[3]=0;
		disk_is_checkpoint[0]=1;
		disk_is_checkpoint[1]=1;
		disk_is_checkpoint[2]=0;
		disk_is_checkpoint[3]=0;
		break;
	case 16:
		*nb_expected_vars=3;
		strncpy(expected_vars[0], "Velocity:9:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[1], "Velocity:15:0:Forward", ADJ_NAME_LEN);
		strncpy(expected_vars[2], "Velocity:16:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_has_value[1]=0;
		memory_has_value[2]=1;
		memory_is_checkpoint[0]=0;
		memory_is_checkpoint[1]=0;
		memory_is_checkpoint[2]=1;
		disk_has_value[0]=1;
		disk_has_value[1]=1;
		disk_has_value[2]=0;
		disk_is_checkpoint[0]=1;
		disk_is_checkpoint[1]=1;
		disk_is_checkpoint[2]=0;
		break;
	case 15:
		*nb_expected_vars=5;
		strncpy(expected_vars[0], "Velocity:9:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_is_checkpoint[0]=0;
		disk_has_value[0]=1;
		disk_is_checkpoint[0]=1;

		strncpy(expected_vars[1], "Velocity:12:0:Forward", ADJ_NAME_LEN);
		memory_has_value[1]=0;
		memory_is_checkpoint[1]=0;
		disk_has_value[1]=1;
		disk_is_checkpoint[1]=1;

		strncpy(expected_vars[2], "Velocity:12:0:Forward", ADJ_NAME_LEN);
		memory_has_value[2]=1;
		memory_is_checkpoint[2]=1;
		disk_has_value[2]=0;
		disk_is_checkpoint[2]=0;

		strncpy(expected_vars[3], "Velocity:14:0:Forward", ADJ_NAME_LEN);
		memory_has_value[3]=1;
		memory_is_checkpoint[3]=1;
		disk_has_value[3]=0;
		disk_is_checkpoint[3]=0;

		strncpy(expected_vars[4], "Velocity:15:0:Forward", ADJ_NAME_LEN);
		memory_has_value[4]=0;
		memory_is_checkpoint[4]=0;
		disk_has_value[4]=1;
		disk_is_checkpoint[4]=1;
		break;

	case 14:
		*nb_expected_vars=5;
		strncpy(expected_vars[0], "Velocity:9:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_is_checkpoint[0]=0;
		disk_has_value[0]=1;
		disk_is_checkpoint[0]=1;

		strncpy(expected_vars[1], "Velocity:12:0:Forward", ADJ_NAME_LEN);
		memory_has_value[1]=0;
		memory_is_checkpoint[1]=0;
		disk_has_value[1]=1;
		disk_is_checkpoint[1]=1;

		strncpy(expected_vars[2], "Velocity:12:0:Forward", ADJ_NAME_LEN);
		memory_has_value[2]=1;
		memory_is_checkpoint[2]=1;
		disk_has_value[2]=0;
		disk_is_checkpoint[2]=0;

		strncpy(expected_vars[3], "Velocity:13:0:Forward", ADJ_NAME_LEN);
		memory_has_value[3]=1;
		memory_is_checkpoint[3]=1;
		disk_has_value[3]=0;
		disk_is_checkpoint[3]=0;

		strncpy(expected_vars[4], "Velocity:14:0:Forward", ADJ_NAME_LEN);
		memory_has_value[4]=1;
		memory_is_checkpoint[4]=1;
		disk_has_value[4]=0;
		disk_is_checkpoint[4]=0;
		break;

	case 5:
		*nb_expected_vars=5;
		strncpy(expected_vars[0], "Velocity:3:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_is_checkpoint[0]=0;
		disk_has_value[0]=1;
		disk_is_checkpoint[0]=1;

		strncpy(expected_vars[1], "Velocity:3:0:Forward", ADJ_NAME_LEN);
		memory_has_value[1]=1;
		memory_is_checkpoint[1]=1;
		disk_has_value[1]=0;
		disk_is_checkpoint[1]=0;

		strncpy(expected_vars[2], "Velocity:4:0:Forward", ADJ_NAME_LEN);
		memory_has_value[2]=1;
		memory_is_checkpoint[2]=1;
		disk_has_value[2]=0;
		disk_is_checkpoint[2]=0;

		strncpy(expected_vars[3], "Velocity:4:0:Forward", ADJ_NAME_LEN);
		memory_has_value[3]=0;
		memory_is_checkpoint[3]=0;
		disk_has_value[3]=1;
		disk_is_checkpoint[3]=1;

		strncpy(expected_vars[4], "Velocity:5:0:Forward", ADJ_NAME_LEN);
		memory_has_value[4]=1;
		memory_is_checkpoint[4]=1;
		disk_has_value[4]=0;
		disk_is_checkpoint[4]=0;
		break;

	case 4:
		*nb_expected_vars=4;
		strncpy(expected_vars[0], "Velocity:3:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_is_checkpoint[0]=0;
		disk_has_value[0]=1;
		disk_is_checkpoint[0]=1;

		strncpy(expected_vars[1], "Velocity:3:0:Forward", ADJ_NAME_LEN);
		memory_has_value[1]=1;
		memory_is_checkpoint[1]=1;
		disk_has_value[1]=0;
		disk_is_checkpoint[1]=0;

		strncpy(expected_vars[2], "Velocity:4:0:Forward", ADJ_NAME_LEN);
		memory_has_value[2]=1;
		memory_is_checkpoint[2]=1;
		disk_has_value[2]=0;
		disk_is_checkpoint[2]=0;

		strncpy(expected_vars[3], "Velocity:4:0:Forward", ADJ_NAME_LEN);
		memory_has_value[3]=0;
		memory_is_checkpoint[3]=0;
		disk_has_value[3]=1;
		disk_is_checkpoint[3]=1;

		break;

	case 3:
		*nb_expected_vars=7;
		strncpy(expected_vars[0], "Velocity:0:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_is_checkpoint[0]=0;
		disk_has_value[0]=1;
		disk_is_checkpoint[0]=1;

		strncpy(expected_vars[1], "Velocity:0:0:Forward", ADJ_NAME_LEN);
		memory_has_value[1]=1;
		memory_is_checkpoint[1]=1;
		disk_has_value[1]=0;
		disk_is_checkpoint[1]=0;

		strncpy(expected_vars[2], "Velocity:1:0:Forward", ADJ_NAME_LEN);
		memory_has_value[2]=0;
		memory_is_checkpoint[2]=0;
		disk_has_value[2]=1;
		disk_is_checkpoint[2]=1;

		strncpy(expected_vars[3], "Velocity:1:0:Forward", ADJ_NAME_LEN);
		memory_has_value[3]=1;
		memory_is_checkpoint[3]=1;
		disk_has_value[3]=0;
		disk_is_checkpoint[3]=0;

		strncpy(expected_vars[4], "Velocity:2:0:Forward", ADJ_NAME_LEN);
		memory_has_value[4]=1;
		memory_is_checkpoint[4]=1;
		disk_has_value[4]=0;
		disk_is_checkpoint[4]=0;

		strncpy(expected_vars[5], "Velocity:3:0:Forward", ADJ_NAME_LEN);
		memory_has_value[5]=0;
		memory_is_checkpoint[5]=0;
		disk_has_value[5]=1;
		disk_is_checkpoint[5]=1;

		strncpy(expected_vars[6], "Velocity:3:0:Forward", ADJ_NAME_LEN);
		memory_has_value[6]=1;
		memory_is_checkpoint[6]=1;
		disk_has_value[6]=0;
		disk_is_checkpoint[6]=0;

		break;

	case 2:
		*nb_expected_vars=5;
		strncpy(expected_vars[0], "Velocity:0:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_is_checkpoint[0]=0;
		disk_has_value[0]=1;
		disk_is_checkpoint[0]=1;

		strncpy(expected_vars[1], "Velocity:0:0:Forward", ADJ_NAME_LEN);
		memory_has_value[1]=1;
		memory_is_checkpoint[1]=1;
		disk_has_value[1]=0;
		disk_is_checkpoint[1]=0;

		strncpy(expected_vars[2], "Velocity:1:0:Forward", ADJ_NAME_LEN);
		memory_has_value[2]=0;
		memory_is_checkpoint[2]=0;
		disk_has_value[2]=1;
		disk_is_checkpoint[2]=1;

		strncpy(expected_vars[3], "Velocity:1:0:Forward", ADJ_NAME_LEN);
		memory_has_value[3]=1;
		memory_is_checkpoint[3]=1;
		disk_has_value[3]=0;
		disk_is_checkpoint[3]=0;

		strncpy(expected_vars[4], "Velocity:2:0:Forward", ADJ_NAME_LEN);
		memory_has_value[4]=1;
		memory_is_checkpoint[4]=1;
		disk_has_value[4]=0;
		disk_is_checkpoint[4]=0;

		break;

	case 1:
		*nb_expected_vars=4;
		strncpy(expected_vars[0], "Velocity:0:0:Forward", ADJ_NAME_LEN);
		memory_has_value[0]=0;
		memory_is_checkpoint[0]=0;
		disk_has_value[0]=1;
		disk_is_checkpoint[0]=1;

		strncpy(expected_vars[1], "Velocity:0:0:Forward", ADJ_NAME_LEN);
		memory_has_value[1]=1;
		memory_is_checkpoint[1]=1;
		disk_has_value[1]=0;
		disk_is_checkpoint[1]=0;

		strncpy(expected_vars[2], "Velocity:1:0:Forward", ADJ_NAME_LEN);
		memory_has_value[2]=0;
		memory_is_checkpoint[2]=0;
		disk_has_value[2]=1;
		disk_is_checkpoint[2]=1;

		strncpy(expected_vars[3], "Velocity:1:0:Forward", ADJ_NAME_LEN);
		memory_has_value[3]=1;
		memory_is_checkpoint[3]=1;
		disk_has_value[3]=0;
		disk_is_checkpoint[3]=0;

		break;

	case 0:
		*nb_expected_vars=0;
		break;

	default:
		printf("Unknown timestep in get_expected_values()");
	  break;
	}
}

/* Check the checkpoint variables */
int test_checkpoints(adj_adjointer *adjointer, int nb_expected_vars, char expected_vars[][ADJ_NAME_LEN], int* memory_has_value, int* memory_is_checkpoint, int* disk_has_value, int* disk_is_checkpoint)
{
	int ierr, i, nb_matched_variables, found_match;
	adj_variable var;
  adj_variable_data* data_ptr;
  char var_name[ADJ_NAME_LEN];

	/* We check that all recorded variables occur in the expectation list.
	 * For that, we loop over all recorded variables.
	 * If the variable name matches with a variable name in the expectation list,
	 * we check if the the storage and checkpoint type match.
	 */
  data_ptr=adjointer->vardata.firstnode;
  while (data_ptr!=NULL)
  {
  	if (data_ptr->equation<0) /* Ignore adjoint variables */
  	{
  		data_ptr=data_ptr->next;
  		continue;
  	}
  	if ((data_ptr->storage.storage_memory_has_value!=ADJ_TRUE) &&
  		  (data_ptr->storage.storage_disk_has_value!=ADJ_TRUE)) /* Ignore variables that are not recorded */
  	{
  		data_ptr=data_ptr->next;
  		continue;
  	}

  	var = adjointer->equations[data_ptr->equation].variable;
  	ierr = adj_variable_str(var, var_name, ADJ_NAME_LEN);
  	if (ierr != ADJ_OK) return ierr;

  	found_match=0;

  	for (i=0; i<nb_expected_vars; i++)
  	{
			if (strcmp(expected_vars[i], var_name)==0)
			{
				if ((memory_has_value[i]==1) && (data_ptr->storage.storage_memory_has_value==ADJ_TRUE))
				{
					if (((memory_is_checkpoint[i]==1) && (data_ptr->storage.storage_memory_is_checkpoint==ADJ_TRUE)) ||
							((memory_is_checkpoint[i]==0) && (data_ptr->storage.storage_memory_is_checkpoint==ADJ_FALSE)))
					{
						found_match=1;
						break;
					}
				}

				if ((disk_has_value[i]==1) && (data_ptr->storage.storage_disk_has_value==ADJ_TRUE))
				{
					if (((disk_is_checkpoint[i]==1) && (data_ptr->storage.storage_disk_is_checkpoint==ADJ_TRUE)) ||
							((disk_is_checkpoint[i]==0) && (data_ptr->storage.storage_disk_is_checkpoint==ADJ_FALSE)))
					{
						found_match=1;
						break;
					}
				}

			}
  	}

  	/* We expect to find all recorded variables in the expectation list */
		if (found_match==0)
		{
			printf("Found an unexpected recorded variable %s\n", var_name);
			adj_test_assert(found_match==1, "Should have worked");
			return ADJ_ERR_NEED_VALUE;
		}

		data_ptr=data_ptr->next;
  }

	/* At this point we know that the recorded variables is a subset of
	 * the expected variables. Next we show that the number of recorded and
	 * expected variables equal, from which we now that both sets are
	 * equivalent.
	 */
	data_ptr=adjointer->vardata.firstnode;
  nb_matched_variables=0; /* Number of checkpoint variables in adjointer */

  while (data_ptr!=NULL)
  {
  	if (data_ptr->equation<0) /* Ignore adjoint variables */
  	{
  		data_ptr=data_ptr->next;
  		continue;
  	}

  	if (data_ptr->storage.storage_memory_has_value)
  		nb_matched_variables++;
  	if (data_ptr->storage.storage_disk_has_value)
  	  		nb_matched_variables++;

  	data_ptr=data_ptr->next;
  }

	/* We expect to the number of expected and recorded variables to be the same */
	if (nb_expected_vars!=nb_matched_variables)
	{
		printf("I expected %i checkpoint variables, but found %i", nb_expected_vars, nb_matched_variables);
		adj_test_assert(nb_expected_vars==nb_matched_variables, "Should have worked");
		return ADJ_ERR_NEED_VALUE;
	}

	return ADJ_OK;
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

void drag_derivative(adj_adjointer* adjointer, adj_variable derivative, int ndepends, adj_variable* variables, adj_vector* dependencies, char* name, adj_vector* output)
{
  int dim=2;
  Vec *output_vec;
  (void) ndepends;
  (void) variables;
  (void) dependencies;
  (void) adjointer;
  (void) derivative;
  (void) name;

  output_vec = (Vec*) malloc(sizeof(Vec));
  VecCreateSeq(PETSC_COMM_SELF, dim, output_vec);
  VecSet(*output_vec, 1.0);
  *output = petsc_vec_to_adj_vector(output_vec);
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
