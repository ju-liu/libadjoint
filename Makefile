###############################################################################
# Find all the relevant objects                                               #
###############################################################################
FSRC = $(wildcard src/*.F90)
FOBJ = $(patsubst src/%.F90,obj/%.o,$(FSRC))

CSRC = $(wildcard src/*.c)
COBJ = $(patsubst src/%.c,obj/%.o,$(CSRC))

###############################################################################
# Variables for unit tests                                                    #
###############################################################################
DISABLED_TESTS = 
FTEST = $(filter-out $(DISABLED_TESTS), $(patsubst src/tests/%,bin/tests/%,$(basename $(filter-out src/tests/test_main.F90, $(wildcard src/tests/*.F90)))))
CTEST = $(filter-out $(DISABLED_TESTS), $(patsubst src/tests/%,bin/tests/%,$(basename $(filter-out src/tests/test_main.c, $(wildcard src/tests/*.c)))))

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
	COMPILER_CFLAGS := -Wall -Wextra -Wunused-parameter -Wunused-result -Wunsafe-loop-optimizations -Wpointer-arith -Wstrict-prototypes -ggdb3 -fstack-protector-all
endif
ifneq (,$(findstring icc, $(CC_VERSION)))
	# icc-specific settings here
	COMPILER_CFLAGS := -Wall 
endif

CFLAGS := $(CFLAGS) $(DBGFLAGS) $(PICFLAG) $(PETSC_CPPFLAGS) -Iinclude/ $(COMPILER_CFLAGS)

###############################################################################
# CPP compiler variables                                                        #
###############################################################################
# Identify CPP compiler
ifeq ($(origin CPP),default)
	CC := mpic++
endif

# Compiler-specific stuff here
CPP_VERSION = $(shell $(CPP) --version 2>&1) $(shell $(CPP) -V 2>&1)
ifneq (,$(findstring g++, $(CPP_VERSION)))
	# g++-specific settings here
	COMPILER_CPPFLAGS := -Wall -Wextra -Wunused-parameter -Wunused-result -Wunsafe-loop-optimizations -Wpointer-arith -Wstrict-prototypes -ggdb3 -fstack-protector-all
endif
ifneq (,$(findstring icc, $(CC_VERSION)))
	# i++-specific settings here
	COMPILER_CPPFLAGS := -Wall 
endif

CPPFLAGS := $(CPPFLAGS) $(DBGFLAGS) $(PICFLAG) $(PETSC_CPPFLAGS) -Iinclude/ $(COMPILER_CFLAGS)

###############################################################################
# F90 compiler variables                                                      #
###############################################################################
# Identify Fortran compiler
ifeq ($(origin FC),default)
	FC := mpif90
endif

# Compiler-specific stuff here
FC_VERSION = $(shell $(FC) --version 2>&1) $(shell $(FC) -V 2>&1)
ifneq (,$(findstring GNU Fortran, $(FC_VERSION)))
	# gfortran-specific settings here
	ifneq (,$(findstring 4.6, $(FC_VERSION)))
		COMPILER_FFLAGS = -fno-whole-file
	endif
endif

ifneq (,$(findstring NAG, $(FC_VERSION)))
	# nag-specific settings here
	COMPILER_FFLAGS = -f2003
endif

FFLAGS := $(FFLAGS) $(DBGFLAGS) $(PICFLAG) $(PETSC_CPPFLAGS) -Iinclude/ -Iinclude/libadjoint $(COMPILER_FFLAGS)

###############################################################################
# Variables for make install                                                  #
###############################################################################
ifeq ($(origin DESTDIR),undefined)
	DESTDIR := /
endif

ifeq ($(origin prefix),undefined)
	prefix := usr/local
endif

###############################################################################
# Variables for static and shared libraries                                   #
###############################################################################
AR = ar
ARFLAGS = cr

LD := $(FC)
LDFLAGS := -shared -Wl,-soname,libadjoint.so

###############################################################################
# Variables for the python bindings                                           #
###############################################################################
H2XML = $(shell which h2xml 2>/dev/null)
XML2PY = $(shell which xml2py 2>/dev/null)
PYDIR = $(shell python -c  "import distutils.sysconfig; print distutils.sysconfig.get_python_lib().replace('/usr/', '$(DESTDIR)/$(prefix)/')")

###############################################################################
# The targets                                                                 #
###############################################################################
all: revolve lib/libadjoint.a lib/libadjoint.so

revolve:
	@echo "  Compiling revolve"
	@cd src/revolve/; make

bin/tests/%: src/tests/%.c src/tests/test_main.c lib/libadjoint.a
	@echo "  CC $@"
	@$(CC) $(CFLAGS) -DTESTNAME=$(notdir $@) -o $@ $< src/tests/test_main.c lib/libadjoint.a $(PETSC_LDFLAGS) $(LIBS)

bin/tests/%: src/tests/%.F90 src/tests/test_main.F90 lib/libadjoint.a
	@echo "  FC $@"
	@$(FC) $(FFLAGS) -DTESTNAME=$(notdir $@) -o $@ $< src/tests/test_main.F90 lib/libadjoint.a $(PETSC_LDFLAGS) $(LIBS)

obj/%.o: src/%.F90
	@echo "  FC $<"
	@$(FC) $(FFLAGS) -c -o $@ $<
	@mv *.mod include/libadjoint 2>/dev/null || true

obj/%.o: src/%.c
	@echo "  CC $<"
	@$(CC) $(CFLAGS) -c -o $@ $<

obj/%.o: src/%.cpp
	@echo "  CPP $<"
	@$(CPP) $(CPPFLAGS) -c -o $@ $<

lib/libadjoint.a: $(COBJ) $(FOBJ) $(CPPOBJ)
	@echo "  AR $@"
	@$(AR) $(ARFLAGS) $@ $(FOBJ) $(COBJ) $(CPPOBJ)

lib/libadjoint.so: $(COBJ) $(FOBJ) $(CPPOBJ)
	@echo "  LD $@"
	@$(LD) $(LDFLAGS) -o $@ $(FOBJ) $(COBJ) $(CPPOBJ) $(PETSC_LDFLAGS) $(LIBS)

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
	@rm -f lib/*.so
	@echo "  RM lib/*.so"
	@rm -f lib/*.py
	@echo "  RM lib/*.py"
	@rm -f tags
	@rm -f include/libadjoint/adj_constants_f.h include/libadjoint/adj_error_handling_f.h
	@cd src/revolve; make clean

test: $(FTEST) $(CTEST)
	@echo "  TEST bin/tests"
	@bin/unittest bin/tests

doc: doc/design/design.pdf doc/manual/manual.pdf

doc/design/design.pdf: doc/design/design.tex doc/design/literature.bib
	@echo "  PDF $@"
	@(cd doc/design && pdflatex -interaction batchmode -shell-escape $(notdir $<) && pdflatex -interaction batchmode -shell-escape $(notdir $<) && bibtex design && pdflatex -interaction batchmode -shell-escape $(notdir $<) && pdflatex -interaction batchmode -shell-escape $(notdir $<)) > /dev/null || \
		echo "    pdflatex failed. Maybe you need to install python-pygments?"

doc/manual/manual.pdf: doc/manual/manual.tex doc/manual/literature.bib
	@echo "  PDF $@"
	@(cd doc/manual && pdflatex -interaction batchmode -shell-escape $(notdir $<) && pdflatex -interaction batchmode -shell-escape $(notdir $<) && bibtex manual && pdflatex -interaction batchmode -shell-escape $(notdir $<) && pdflatex -interaction batchmode -shell-escape $(notdir $<)) > /dev/null || \
		echo "    pdflatex failed. Maybe you need to install python-pygments?"

ifneq (,$(CTAGS))
all: tags
test: tags
tags: $(FSRC) $(CSRC)
	@echo "  CTAGS src/*.c src/*.F90"
	@$(CTAGS) src/*.c src/*.F90
endif

ifneq (,$(H2XML))
all: lib/libadjoint.py
test: lib/libadjoint.py
install: lib/libadjoint.py
lib/libadjoint.py: lib/libadjoint.so
	@echo "  H2XML  include/libadjoint/libadjoint.h"
	@$(H2XML) -q -I. include/libadjoint/libadjoint.h -o lib/libadjoint.xml
	@echo "  XML2PY lib/libadjoint.py"
	@$(XML2PY) -r '^adj.*' -l lib/libadjoint.so lib/libadjoint.xml -o lib/libadjoint.py
	@rm -f lib/libadjoint.xml
	@chmod a-x lib/libadjoint.py
endif

install: lib/libadjoint.a lib/libadjoint.so
	@echo "  INSTALL $(DESTDIR)/$(prefix)/lib"
	@install -d $(DESTDIR)/$(prefix)/lib
	@install lib/libadjoint.a $(DESTDIR)/$(prefix)/lib
	@install lib/libadjoint.so $(DESTDIR)/$(prefix)/lib
ifneq (,$(H2XML))
	@echo "  INSTALL $(PYDIR)"
	@install -d $(PYDIR)
	@install -m 644 lib/libadjoint.py $(PYDIR)
	@sed -i "s@CDLL('lib/libadjoint.so')@CDLL('/$(prefix)/lib/libadjoint.so')@" $(PYDIR)/libadjoint.py
endif
	@echo "  INSTALL $(DESTDIR)/$(prefix)/include/libadjoint"
	@install -d $(DESTDIR)/$(prefix)/include/libadjoint
	@install include/libadjoint/* $(DESTDIR)/$(prefix)/include/libadjoint

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

