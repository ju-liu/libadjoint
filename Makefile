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

# Identify if PETSc is installed
PETSC_CPPFLAGS = $(shell make -f cfg/petsc_makefile getincludedirs 2>/dev/null)
PETSC_LDFLAGS  = $(shell make -f cfg/petsc_makefile getlinklibs 2>/dev/null)
ifeq (,$(PETSC_CPPFLAGS))
	PETSC_CPPFLAGS := -UHAVE_PETSC
else
	PETSC_CPPFLAGS := $(PETSC_CPPFLAGS) -DHAVE_PETSC
endif

DBGFLAGS = -g -O0
CFLAGS = $(DBGFLAGS) $(PETSC_CPPFLAGS) -Iinclude/ $(COMPILER_CFLAGS)

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

FFLAGS = $(DBGFLAGS) $(PETSC_CPPFLAGS) -Iinclude/ $(COMPILER_FFLAGS)

AR = ar
ARFLAGS = cr

FSRC = $(wildcard src/*.F90)
FOBJ = $(patsubst src/%.F90,obj/%.o,$(FSRC))

CSRC = $(wildcard src/*.c)
COBJ = $(patsubst src/%.c,obj/%.o,$(CSRC))

DISABLED_TESTS =  
FTEST = $(filter-out $(DISABLED_TESTS), $(patsubst src/tests/%,bin/tests/%,$(basename $(filter-out src/tests/test_main.F90, $(wildcard src/tests/*.F90)))))
CTEST = $(filter-out $(DISABLED_TESTS), $(patsubst src/tests/%,bin/tests/%,$(basename $(filter-out src/tests/test_main.c, $(wildcard src/tests/*.c)))))

CTAGS = $(shell which ctags)

all: lib/libadjoint.a

obj/tests/%: src/tests/%.F90 src/tests/test_main.F90 lib/libadjoint.a
	@echo "  FC $@"
	@$(FC) $(FFLAGS) -DTESTNAME=$(notdir $@) -o $@ $< src/tests/test_main.F90 -Llib/ -ladjoint $(PETSC_LDFLAGS)

bin/tests/%: src/tests/%.c src/tests/test_main.c lib/libadjoint.a
	@echo "  CC $@"
	@$(CC) $(CFLAGS) -DTESTNAME=$(notdir $@) -o $@ $< src/tests/test_main.c -Llib/ -ladjoint -I./include/libadjoint $(PETSC_LDFLAGS)

obj/%.o: src/%.F90
	@echo "  FC $<"
	@$(FC) $(FFLAGS) -c -o $@ $<
	@mv *.mod include/libadjoint 2>/dev/null || true

obj/%.o: src/%.c
	@echo "  CC $<"
	@$(CC) $(CFLAGS) -c -o $@ $<

lib/libadjoint.a: $(COBJ) $(FOBJ)
	@echo "  AR $@"
	@$(AR) $(ARFLAGS) $@ $(FOBJ) $(COBJ)

clean:
	@rm -f obj/*.o
	@echo "  RM obj/*.o"
	@rm -f obj/tests/*.o
	@echo "  RM obj/tests/*.o"
	@rm -f bin/tests/*
	@echo "  RM bin/tests/*"
	@rm -f include/libadjoint/*.mod include/*.mod *.mod
	@echo "  RM include/libadjoint/*.mod"
	@rm -f doc/*.pdf doc/*.bbl doc/*.log doc/*.aux doc/*.blg doc/*.spl doc/*.brf doc/*.out doc/*.toc
	@echo "  RM doc/*.pdf"
	@rm -f lib/*.a
	@echo "  RM lib/*.a"
	@rm -f tags
	@rm -f include/libadjoint/adj_constants_f.h

test: $(FTEST) $(CTEST)
	@echo "  TEST bin/tests"
	@bin/unittest bin/tests

doc: doc/design.pdf

doc/design.pdf: doc/design.tex doc/literature.bib
	@echo "  PDF $@"
	@(cd doc && pdflatex -interaction batchmode -shell-escape $(notdir $<) && pdflatex -interaction batchmode -shell-escape $(notdir $<) && bibtex design && pdflatex -interaction batchmode -shell-escape $(notdir $<) && pdflatex -interaction batchmode -shell-escape $(notdir $<)) > /dev/null || \
		echo "    pdflatex failed. Maybe you need to install python-pygments?"

ifneq (,$(CTAGS))
lib/libadjoint.a: tags
tags: $(FSRC) $(CSRC)
	@echo "  CTAGS src/*.c src/*.F90"
	@$(CTAGS) src/*.c src/*.F90
endif

include/libadjoint/adj_fortran.h: include/libadjoint/adj_constants_f.h include/libadjoint/adj_error_handling_f.h
# replace C comments with F90 comments
include/libadjoint/adj_constants_f.h: include/libadjoint/adj_constants.h
	@echo "  SED $@"
	@sed -e 's@/\*@!@' -e 's@\*/@@' -e 's@ADJ_CONSTANTS_H@ADJ_CONSTANTS_F_H@' $< > $@
include/libadjoint/adj_error_handling_f.h: include/libadjoint/adj_error_handling.h
	@echo "  SED $@"
	@grep '^#' $< | grep -v '^#include' > $@

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

