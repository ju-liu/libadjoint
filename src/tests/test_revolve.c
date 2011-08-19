#include <stdlib.h>
#include <stdio.h>
#include "libadjoint/revolve_c.h"
#include "libadjoint/adj_test_main.h"
#include "libadjoint/adj_test_tools.h"

void test_revolve(void)
{
  int steps=5;
  int snaps;
  int snaps_in_ram;
  int timestep;

  int checkpoint_time[2]; /* Saves the timestep of the checkpoints */

  CACTION action;
  CRevolve r;

  snaps = revolve_adjust(r, steps);
  printf("Doing %i timesteps with %i checkpoints\n", steps, snaps);
  adj_test_assert(snaps == 2, "Snaps should be 2");

  r = revolve_create_offline(steps, snaps);

  timestep = 0;
  do
  {
    action = revolve(r);
    printf("Action: %s\n", revolve_caction_string(action));
    switch(action)
    {
      case CACTION_TAKESHOT:
        printf("Save a checkpoint in slot %i\n", revolve_getcheck(r)); 
        printf("------------------------------\n");
        checkpoint_time[revolve_getcheck(r)] = timestep;
        break;
      case CACTION_ADVANCE:
        printf("Advance from timestep %i to %i\n", revolve_getoldcapo(r), revolve_getcapo(r));
        printf("------------------------------\n");
        timestep = revolve_getcapo(r);
        break;
      case CACTION_FIRSTRUN:
        printf("Advance from timestep %i to %i\n", timestep, timestep+1); 
        printf("Compute adjoint from timestep %i to %i\n", timestep+1, timestep);
        printf("------------------------------\n");
        break;
      case CACTION_YOUTURN:
        printf("Advance from timestep %i to %i\n", timestep, timestep+1); 
        printf("Compute adjoint from timestep %i to %i\n", timestep+1, timestep);
        printf("------------------------------\n");
        timestep = timestep - 1;
        break;
      case CACTION_RESTORE:
        printf("Restore checkpoint in slot %i containing timestep %i\n", revolve_getcheck(r), checkpoint_time[revolve_getcheck(r)]); 
        printf("------------------------------\n");
        break;
      case CACTION_ERROR:
        adj_test_assert(1==0, "Irregular termination of revolve");
        break;
      case CACTION_TERMINATE:
        break;
      default:
        adj_test_assert(1==0, "Unknown revolve action");
        break;
    }
  }
  while ((action!=CACTION_TERMINATE) && (action!=CACTION_ERROR));

  revolve_destroy(r);

}
