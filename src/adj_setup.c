#ifdef HAVE_PETSC
#include "petscsys.h"  
#include "petsc.h"
#endif
  
void  adj_init(int ierr)
{
  ierr=0;
#ifdef HAVE_PETSC
  int argc=0;
  PetscInitialize(&argc, PETSC_NULL, PETSC_NULL, PETSC_NULL);
#endif
}

void adj_finalize(int ierr)
{
  ierr=0;
#ifdef HAVE_PETSC
  PetscFinalize();
#endif
}
