#include "petscversion.h"
#ifdef HAVE_PETSC_MODULES
#if PETSC_VERSION_MINOR==0
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscpcdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscisdef.h"
#else
#include "finclude/petscdef.h"
#endif
#else
#include "finclude/petsc.h"
#if PETSC_VERSION_MINOR==0
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscviewer.h"
#include "finclude/petscis.h"
#endif
#endif
