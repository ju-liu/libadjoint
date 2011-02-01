#include "petscversion.h"
#ifdef HAVE_PETSC_MODULES
#if PETSC_VERSION_MINOR==0
#include "petscvecdef.h"
#include "petscmatdef.h"
#include "petsckspdef.h"
#include "petscpcdef.h"
#include "petscviewerdef.h"
#include "petscisdef.h"
#else
#include "petscdef.h"
#endif
#else
#include "petsc.h"
#if PETSC_VERSION_MINOR==0
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscpc.h"
#include "petscviewer.h"
#include "petscis.h"
#endif
#endif

