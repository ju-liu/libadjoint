###############################################################################
# Print out a helpful table of dependencies                                   #
###############################################################################
DEPENDENCIES = $(shell bin/checkdeps)

###############################################################################
# Find all the relevant objects                                               #
###############################################################################
FSRC = $(wildcard src/*.F90)
FOBJ = $(patsubst src/%.F90,obj/%.o,$(FSRC))

CSRC = $(wildcard src/*.c)
COBJ = $(patsubst src/%.c,obj/%.o,$(CSRC))

CXXSRC = $(wildcard src/*.cpp)
CXXOBJ = $(patsubst src/%.cpp,obj/%.o,$(CXXSRC))

###############################################################################
# Variables for unit tests                                                    #
###############################################################################
DISABLED_TESTS = 
FTEST = $(filter-out $(DISABLED_TESTS), $(patsubst src/tests/%,bin/tests/%,$(basename $(filter-out src/tests/test_main.F90, $(wildcard src/tests/*.F90)))))
CTEST = $(filter-out $(DISABLED_TESTS), $(patsubst src/tests/%,bin/tests/%,$(basename $(filter-out src/tests/test_main.c, $(wildcard src/tests/*.c)))))
PTEST = $(filter-out $(DISABLED_TESTS), $(patsubst src/tests/%,bin/tests/%,$(basename $(wildcard src/tests/*.py))))

###############################################################################
# Variables for unit tests                                                    #
###############################################################################
CTAGS = $(shell which ctags 2>/dev/null)

###############################################################################
# Generic compiler flags                                                      #
###############################################################################
# Identify if PETSc is installed
MAKE = make -s
PETSC_CPPFLAGS = $(shell $(MAKE) -f cfg/petsc_makefile getincludedirs 2>/dev/null)
PETSC_LDFLAGS  = $(shell $(MAKE) -f cfg/petsc_makefile getlinklibs 2>/dev/null)
ifeq (,$(PETSC_CPPFLAGS))
	PETSC_CPPFLAGS := # want to have -UHAVE_PETSC, but that causes confusion on some fortran compilers (e.g. nag) and it isn't really necessary
else
	PETSC_CPPFLAGS := $(PETSC_CPPFLAGS) -DHAVE_PETSC
endif

# Identify if SLEPc is installed
ifneq (,$(SLEPC_DIR))
	SLEPC_CPPFLAGS := -I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include -DHAVE_SLEPC
	SLEPC_LDFLAGS  := -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -lslepc
else
	SLEPC_CPPFLAGS :=
	SLEPC_LDFLAGS  :=
endif

DBGFLAGS = -g -O0
PICFLAG = -fPIC

###############################################################################
# C compiler variables                                                        #
###############################################################################
# Identify C compiler
ifeq ($(origin CC),default)
	CC := mpicc
endif

# Compiler-specific stuff here
CC_VERSION = $(shell $(CC) --version 2>&1) $(shell $(CC) -V 2>&1)
ifneq (,$(findstring gcc, $(CC_VERSION)))
	# gcc-specific settings here
	COMPILER_CFLAGS := -Wall -Wextra -Wunused-parameter -Wunsafe-loop-optimizations -Wpointer-arith -Wstrict-prototypes -ggdb3 -fstack-protector-all -lstdc++
endif
ifneq (,$(findstring icc, $(CC_VERSION)))
	# icc-specific settings here
	COMPILER_CFLAGS := -Wall 
endif

CFLAGS := $(CFLAGS) $(DBGFLAGS) $(PICFLAG) $(SLEPC_CPPFLAGS) $(PETSC_CPPFLAGS) -Iinclude/ $(COMPILER_CFLAGS)

###############################################################################
# CXX compiler variables                                                        #
###############################################################################
ifeq ($(origin CXX),default)
	CXX := mpic++
endif

# Compiler-specific stuff Vere
CXX_VERSION = $(shell $(CXX) --version 2>&1) $(shell $(CXX) -V 2>&1)
ifneq (,$(findstring g++, $(CXX_VERSION)))
	# g++-specific settings here
	COMPILER_CXXFLAGS := -Wunsafe-loop-optimizations -Wpointer-arith -ggdb3 -fstack-protector-all
endif
ifneq (,$(findstring icpc, $(CXX_VERSION)))
	# i++-specific settings here
	COMPILER_CXXFLAGS := 
endif

CXXFLAGS := $(CXXFLAGS) $(DBGFLAGS) $(PICFLAG) $(SLEPC_CPPFLAGS) $(PETSC_CPPFLAGS) -Iinclude/ $(COMPILER_CXXFLAGS)

###############################################################################
# F90 compiler variables                                                      #
###############################################################################
# Identify Fortran compiler
ifeq ($(origin FC),default)
  ifneq ($(shell which mpif90),)
	FC := mpif90
  else
    FC :=
  endif
endif

# Compiler-specific stuff here
COMPILER_FFLAGS := 
ifneq (,$(FC))
FC_VERSION := $(shell $(FC) --version 2>&1) $(shell $(FC) -V 2>&1)
else
FC_VERSION := 
endif

ifneq (,$(findstring GNU Fortran, $(FC_VERSION)))
	# gfortran-specific settings here
        COMPILER_FFLAGS := $(COMPILER_FFLAGS) -lstdc++
        ifneq (,$(findstring 4.4, $(FC_VERSION)))
                FC :=
        endif
	ifneq (,$(findstring 4.6, $(FC_VERSION)))
		COMPILER_FFLAGS := $(COMPILER_FFLAGS) -fno-whole-file
	endif
endif

ifneq (,$(findstring NAG, $(FC_VERSION)))
	# nag-specific settings here
	COMPILER_FFLAGS := $(COMPILER_FFLAGS) -f2003
endif

FFLAGS := $(FFLAGS) $(DBGFLAGS) $(PICFLAG) $(SLEPC_CPPFLAGS) $(PETSC_CPPFLAGS) -Iinclude/ -Iinclude/libadjoint $(COMPILER_FFLAGS)

ifneq (,$(FC))
OBJECTS := $(COBJ) $(FOBJ) $(CXXOBJ)
else
OBJECTS := $(COBJ) $(CXXOBJ)
endif

###############################################################################
# Variables for make install                                                  #
###############################################################################
ifeq ($(origin DESTDIR),undefined)
	DESTDIR := /
endif
ABSDESTDIR := $(shell bin/abspath $(DESTDIR))

ifeq ($(origin prefix),undefined)
	prefix := usr/local
endif

###############################################################################
# Variables for static and shared libraries                                   #
###############################################################################
AR = ar
ARFLAGS = cr

ifneq (,$(FC))
LD := $(FC)
else
LD := $(CXX)
endif

ifeq (,$(findstring Darwin, $(shell uname -a)))
SLIB := libadjoint.so
CXXLIBS := -lstdc++ -lsupc++
LDFLAGS := $(CXXLIBS) -shared -Wl,-soname,$(SLIB)
SEDFLG  := 
else
SLIB := libadjoint.dylib
CXXLIBS := -lstdc++ -lsupc++
LDFLAGS := $(CXXLIBS) -shared -Wl,-install_name,$(SLIB)
SEDFLG  := ""
endif

###############################################################################
# Variables for the python bindings                                           #
###############################################################################
ifeq (,$(findstring Darwin, $(shell uname -a)))
CPP := cpp
else
CPP := gcc -E -U__BLOCKS__
endif
GCCXML = $(shell which gccxml)
H2XML = python/ctypeslib/scripts/h2xml.py
XML2PY = python/ctypeslib/scripts/xml2py.py
PYDIR = $(shell python -c  "import distutils.sysconfig; print distutils.sysconfig.get_python_lib().replace('/usr/', '$(ABSDESTDIR)/$(prefix)/')")


###############################################################################
# The targets                                                                 #
###############################################################################
all: lib/libadjoint.a $(SLIB)

bin/tests/%: src/tests/%.c src/tests/test_main.c lib/libadjoint.a
	@echo "  CC $@"
	@$(CC) $(CFLAGS) -DTESTNAME=$(notdir $@) -o $@ $< src/tests/test_main.c lib/libadjoint.a $(SLEPC_LDFLAGS) $(PETSC_LDFLAGS) $(LIBS) $(CXXLIBS)

ifneq ($(FC),)
bin/tests/%: src/tests/%.F90 src/tests/test_main.F90 lib/libadjoint.a
	@echo "  FC $@"
	@$(FC) $(FFLAGS) -DTESTNAME=$(notdir $@) -o $@ $< src/tests/test_main.F90 lib/libadjoint.a $(SLEPC_LDFLAGS) $(PETSC_LDFLAGS) $(LIBS) $(CXXLIBS)
endif

bin/tests/%: src/tests/%.py pybuild
	@echo "  PY $@"
	@cp src/tests/$(notdir $@).py bin/tests

ifneq ($(FC),)
obj/%.o: src/%.F90
	@echo "  FC $<"
	@$(FC) $(FFLAGS) -c -o $@ $<
	@mv *.mod include/libadjoint 2>/dev/null || true
endif

obj/%.o: src/%.c
	@echo "  CC $<"
	@$(CC) $(CFLAGS) -c -o $@ $<

obj/%.o: src/%.cpp
	@echo "  C++ $<"
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

.PHONY: checkdeps
checkdeps:
	@echo $(DEPENDENCIES) > /dev/null

lib/libadjoint.a: $(OBJECTS)
	@echo "  AR $@"
	@$(AR) $(ARFLAGS) $@ obj/*.o

$(SLIB): $(OBJECTS)
	@echo "  LD $@"
	@$(LD) -o $@ obj/*.o $(SLEPC_LDFLAGS) $(PETSC_LDFLAGS) $(LIBS) $(LDFLAGS)

clean:
	@rm -f obj/*.o
	@echo "  RM obj/*.o"
	@rm -f obj/tests/*.o
	@echo "  RM obj/tests/*.o"
	@rm -f bin/tests/*
	@echo "  RM bin/tests/*"
	@rm -f include/libadjoint/*.mod include/*.mod *.mod
	@echo "  RM include/libadjoint/*.mod"
	@rm -f doc/*/*.pdf doc/*/*.bbl doc/*/*.log doc/*/*.blg doc/*/*.spl doc/*/*.brf doc/*/*.out doc/*/*.toc
	@echo "  RM doc/*/*.pdf"
	@rm -f lib/*.a
	@echo "  RM lib/*.a"
	@rm -f lib/*.so lib/*.dylib
	@echo "  RM lib/*.so"
	@rm -f python/libadjoint/clibadjoint.py
	@echo "  RM python/libadjoint/clibadjoint.py"
	@rm -f python/libadjoint/python_utils.so python/libadjoint/python_utils.dylib
	@echo "  RM python/libadjoint/python_utils.so"
	@rm -f python/libadjoint/clibadjoint_constants.py
	@echo "  RM python/libadjoint/clibadjoint_constants.py"
	@rm -rf python/build python/libadjoint/*.pyc
	@echo "  RM python/build"
	@rm -f tags
	@rm -f include/libadjoint/adj_constants_f.h include/libadjoint/adj_error_handling_f.h

ifneq (,$(FC))
compiled_tests: $(FTEST) $(CTEST)
else
compiled_tests: $(CTEST)
endif

ifneq (,$(GCCXML))
test: compiled_tests $(PTEST) pybuild
else
test: compiled_tests
endif
	@export PYTHONPATH=$(PWD)/python
	@echo "  TEST bin/tests"
	@bin/unittest bin/tests

pybuild: python/build

python/build: python/libadjoint/*.c
	@echo "  PYBUILD python/libadjoint"
	@cd python && python setup.py -q build
	@bin/link_python

doc: doc/design/design.pdf doc/manual/manual.pdf

doc/design/design.pdf: doc/design/design.tex doc/design/literature.bib
	@echo "  PDF $@"
	@(cd doc/design && pdflatex -interaction batchmode -shell-escape $(notdir $<) && pdflatex -interaction batchmode -shell-escape $(notdir $<) && bibtex design && pdflatex -interaction batchmode -shell-escape $(notdir $<) && pdflatex -interaction batchmode -shell-escape $(notdir $<)) > /dev/null || \
		echo "    pdflatex failed. Maybe you need to install python-pygments?"

doc/manual/manual.pdf: doc/manual/*.tex doc/manual/literature.bib
	@echo "  PDF $@"
	@(cd doc/manual && make > /dev/null 2>&1) || \
		echo "    pdflatex failed. Maybe you need to install python-pygments?"

ifneq (,$(CTAGS))
all: tags
test: tags
tags: $(FSRC) $(CSRC)
	@echo "  CTAGS src/*.c src/*.F90"
	@$(CTAGS) src/*.c src/*.F90
endif

ifneq (,$(GCCXML))
python: python/libadjoint/clibadjoint.py python/libadjoint/clibadjoint_constants.py
all: python
test: python
install: python

python/libadjoint/clibadjoint.py: $(SLIB)
	@echo "  H2XML  include/libadjoint/libadjoint.h"
	@$(CPP) -DPYTHON_BINDINGS include/libadjoint/libadjoint.h > include/libadjoint/pylibadjoint.h
	@sed -i $(SEDFLG) "s/__builtin___stpncpy_chk/__builtin___strncpy_chk/" include/libadjoint/pylibadjoint.h
	@$(H2XML) -q -I. include/libadjoint/pylibadjoint.h -o python/libadjoint/libadjoint.xml
	@rm -f include/libadjoint/pylibadjoint.h
	@echo "  XML2PY python/libadjoint/clibadjoint.py"
	@$(XML2PY) -r '^adj.*' -l $(shell python bin/realpath $(SLIB)) python/libadjoint/libadjoint.xml -o python/libadjoint/clibadjoint.py
	@rm -f python/libadjoint/libadjoint.xml
	@chmod a-x python/libadjoint/clibadjoint.py
python/libadjoint/clibadjoint_constants.py:
	@echo "  PYTHON tools/create_python_constants.py"
	@python ./tools/create_python_constants.py
endif

install: lib/libadjoint.a $(SLIB)
	@echo "  INSTALL $(ABSDESTDIR)/$(prefix)/lib"
	@install -d $(ABSDESTDIR)/$(prefix)/lib
	@install lib/libadjoint.a $(ABSDESTDIR)/$(prefix)/lib
	@install $(SLIB) $(ABSDESTDIR)/$(prefix)/lib
ifneq (,$(GCCXML))
	@echo "  INSTALL $(PYDIR)"
ifeq ($(LIBADJOINT_BUILDING_DEBIAN),yes)
	@cd python; for PYTHON in $(shell pyversions -r); do echo $$PYTHON; $$PYTHON setup.py install --prefix=$(ABSDESTDIR)/$(prefix) $(LIBADJOINT_PYTHON_INSTALL_ARGS); done
else
	@cd python; python setup.py install --prefix=$(ABSDESTDIR)/$(prefix) $(LIBADJOINT_PYTHON_INSTALL_ARGS)
endif
	@find $(ABSDESTDIR)/$(prefix) -name clibadjoint.py | xargs sed -i $(SEDFLG) -e "s@CDLL('$(SLIB)')@CDLL('/$(prefix)/lib/$(SLIB)')@" -e "s@CDLL('$(shell python bin/realpath $(SLIB))')@CDLL('/$(prefix)/lib/$(SLIB)')@"
endif
	@echo "  INSTALL $(ABSDESTDIR)/$(prefix)/include/libadjoint"
	@install -d $(ABSDESTDIR)/$(prefix)/include/libadjoint
	@install include/libadjoint/* $(ABSDESTDIR)/$(prefix)/include/libadjoint

include/libadjoint/adj_fortran.h: include/libadjoint/adj_constants_f.h include/libadjoint/adj_error_handling_f.h
# replace C comments with F90 comments
include/libadjoint/adj_constants_f.h: include/libadjoint/adj_constants.h
	@echo "  SED $@"
	@sed -e 's@/\*@!@' -e 's@\*/@@' -e 's@ADJ_CONSTANTS_H@ADJ_CONSTANTS_F_H@' -e '/adj_scalar /d' $< > $@
include/libadjoint/adj_error_handling_f.h: include/libadjoint/adj_error_handling.h
	@echo "  SED $@"
	@grep '^#' $< | grep -v '^#include' | grep -v 'CHKMALLOC' > $@

# You can generate some of this with: for i in src/*.c; do gcc -Iinclude/ -MM $i; done | sed -e 's@\\@@' -e 's@^adj@obj/adj@'
obj/adj_data_structures.o: include/libadjoint/adj_data_structures.h include/libadjoint/adj_constants.h include/libadjoint/adj_error_handling.h include/libadjoint/uthash.h
obj/adj_error_handling.o: include/libadjoint/adj_error_handling.h
obj/adj_adjointer_routines.o: include/libadjoint/adj_adjointer_routines.h include/libadjoint/adj_data_structures.h include/libadjoint/adj_constants.h \
	                            include/libadjoint/uthash.h include/libadjoint/adj_variable_lookup.h include/libadjoint/adj_error_handling.h
obj/adj_variable_lookup.o: include/libadjoint/adj_variable_lookup.h include/libadjoint/adj_data_structures.h include/libadjoint/adj_constants.h include/libadjoint/uthash.h  include/libadjoint/adj_error_handling.h
obj/adj_petsc_datastructures.o:  include/libadjoint/adj_petsc_datastructures.h include/libadjoint/adj_data_structures.h
obj/adj_evaluation.o: include/libadjoint/adj_evaluation.h include/libadjoint/adj_data_structures.h include/libadjoint/adj_constants.h include/libadjoint/uthash.h include/libadjoint/adj_error_handling.h \
	                    include/libadjoint/adj_adjointer_routines.h include/libadjoint/adj_data_structures.h include/libadjoint/adj_variable_lookup.h include/libadjoint/adj_error_handling.h
obj/adj_petsc_data_structures.o: include/libadjoint/adj_petsc_data_structures.h include/libadjoint/adj_constants.h include/libadjoint/adj_adjointer_routines.h include/libadjoint/adj_data_structures.h include/libadjoint/uthash.h \
	                               include/libadjoint/adj_variable_lookup.h include/libadjoint/adj_error_handling.h
obj/adj_simplification.o: include/libadjoint/adj_simplification.h include/libadjoint/adj_data_structures.h include/libadjoint/adj_constants.h include/libadjoint/uthash.h include/libadjoint/adj_adjointer_routines.h \
	                        include/libadjoint/adj_variable_lookup.h include/libadjoint/adj_error_handling.h
obj/adj_test_tools.o: include/libadjoint/adj_test_tools.h
obj/adj_fortran.o: include/libadjoint/adj_fortran.h
obj/adj_core.o: include/libadjoint/adj_core.h include/libadjoint/adj_data_structures.h include/libadjoint/adj_constants.h include/libadjoint/uthash.h include/libadjoint/adj_error_handling.h include/libadjoint/adj_variable_lookup.h \
	              include/libadjoint/adj_adjointer_routines.h include/libadjoint/adj_evaluation.h include/libadjoint/adj_data_structures.h include/libadjoint/adj_error_handling.h include/libadjoint/adj_adjointer_routines.h \
	              include/libadjoint/adj_simplification.h
obj/adj_dictionary.o: include/libadjoint/adj_dictionary.h include/libadjoint/adj_error_handling.h include/libadjoint/adj_constants.h

