# Identify C compiler
ifeq (,${CC})
	CC := mpicc
endif

# Identify if PETSc is installed
PETSC_CPPFLAGS = $(shell make -f cfg/petsc_makefile getincludedirs 2>/dev/null)
PETSC_LDFLAGS  = $(shell make -f cfg/petsc_makefile getlinklibs 2>/dev/null)
ifeq (,$(PETSC_CPPFLAGS)
	PETSC_CPPFLAGS := -UHAVE_PETSC
else
	PETSC_CPPFLAGS := $(PETSC_CPPFLAGS) -DHAVE_PETSC
endif

DBGFLAGS = -g -O0
CFLAGS = $(DBGFLAGS) $(PETSC_CPPFLAGS) -Iinclude/

# Identify Fortran compiler
ifeq (,${FC})
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

FFLAGS = $(DBGFLAGS) $(PETSC_CPPFLAGS) -Iinclude/ -DHAVE_PETSC $(COMPILER_FFLAGS)

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

obj/tests/%: src/tests/%.c src/tests/test_main.c lib/libadjoint.a
	@echo "  CC $@"
	@$(CC) $(CFLAGS) -DTESTNAME=$(notdir $@) -o $@ $< src/tests/test_main.c -Llib/ -ladjoint $(PETSC_LDFLAGS)

obj/%.o: src/%.F90
	@echo "  FC $<"
	@$(FC) $(FFLAGS) -c -o $@ $<
	@mv *.mod include/ 2>/dev/null || true

obj/%.o: src/%.c
	@echo "  CC $<"
	@$(CC) $(CFLAGS) -c -o $@ $<

lib/libadjoint.a: $(FOBJ) $(COBJ)
	@echo "  AR $@"
	@$(AR) $(ARFLAGS) $@ $(FOBJ) $(COBJ)

clean:
	@rm -f obj/*.o
	@echo "  RM obj/*.o"
	@rm -f obj/tests/*.o
	@echo "  RM obj/tests/*.o"
	@rm -f bin/tests/*
	@echo "  RM bin/tests/*"
	@rm -f include/*.mod
	@echo "  RM include/*.mod"
	@rm -f doc/*.pdf doc/*.bbl doc/*.log doc/*.aux doc/*.blg doc/*.spl doc/*.brf doc/*.out doc/*.toc
	@echo "  RM doc/*.pdf"
	@rm -f lib/*.a
	@echo "  RM lib/*.a"
	@rm -f tags

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
tags: $(FSRC)
	@echo "  CTAGS src/*.F90"
	@$(CTAGS) src/*.F90
endif
