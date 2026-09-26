# fortran/tests/Makefile
# Build and run Fortran binding tests against libhexpr.
#
# Usage:
#   make [FC=gfortran]   # default
#   make FC=flang
#   make FC=ifx          # source oneapi setvars.sh first if needed
#
# Targets: all (default), check, clean

# GNU Make's built-in default for FC is f77; override it here so that
# "make" uses gfortran by default.  Users can still pass FC=flang etc.
ifeq ($(origin FC), default)
  FC := gfortran
endif

ifeq ($(FC),gfortran)
  FCFLAGS := -Wall -Wextra -O2
else ifeq ($(FC),flang)
  FCFLAGS := -Wall -O2
else ifeq ($(FC),ifx)
  FCFLAGS := -warn all -O2
else
  FCFLAGS := -O2
endif

# libhexpr already contains Wigner symbols (linked with --whole-archive);
# no separate -lwigner needed. libhexpr links only against libm.
#
# Platform link: -L../../c -lhexpr resolves libhexpr.so on Linux and
# libhexpr.dylib on macOS automatically. On macOS bake an rpath so the test
# binaries (in fortran/tests/) find the dylib (in c/) without any env var.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  RPATH := -Wl,-rpath,@loader_path/../../c
else
  RPATH :=
endif
LIBS      := -L../../c -lhexpr $(RPATH) -lm

F90_TESTS := test_smoke test_drt test_expr test_csfset test_pair test_completeness
F77_TESTS := test_f77_smoke test_f77_completeness test_f77_eval
TESTS     := $(F90_TESTS) $(F77_TESTS)

.PHONY: all check clean

all: $(TESTS)

# ---- Module (compiled into tests/ directory) ---------------------
# The .mod file lands here; -I. makes it visible when compiling tests.
#
# One compile produces two files, hexpr_mod.o and hexpr_mod.mod, and GNU
# make has no way to say that.  This rule used to read
#
#     hexpr_mod.mod hexpr_mod.o: ../src/hexpr_mod.f90
#
# which make reads as TWO independent rules, not as one rule with two
# outputs.  The tests depended on hexpr_mod.o alone, so with the .o up to
# date the recipe never ran and a missing or stale hexpr_mod.mod was
# never regenerated -- every test then compiled against whatever .mod
# happened to be on disk.  That is silent when the stale interface still
# matches, and reports a missing procedure when it does not.
#
# Both files are therefore produced through one stamp, and the tests
# depend on the stamp rather than on either file.
hexpr_mod.stamp: ../src/hexpr_mod.f90
	$(FC) $(FCFLAGS) -c $< -o hexpr_mod.o
	@touch $@

# If either artifact was removed by hand while the stamp still looks
# current, drop the stamp so the next build regenerates both.  Evaluated
# when the makefile is read, before any target is considered.
_MOD_GUARD := $(shell test -f hexpr_mod.o -a -f hexpr_mod.mod || rm -f hexpr_mod.stamp)

# ---- F90 tests --------------------------------------------------
test_smoke: test_smoke.f90 hexpr_mod.stamp
	$(FC) $(FCFLAGS) -I. $< hexpr_mod.o $(LIBS) -o $@

test_drt: test_drt.f90 hexpr_mod.stamp
	$(FC) $(FCFLAGS) -I. $< hexpr_mod.o $(LIBS) -o $@

test_expr: test_expr.f90 hexpr_mod.stamp
	$(FC) $(FCFLAGS) -I. $< hexpr_mod.o $(LIBS) -o $@

test_csfset: test_csfset.f90 hexpr_mod.stamp
	$(FC) $(FCFLAGS) -I. $< hexpr_mod.o $(LIBS) -o $@

test_pair: test_pair.f90 hexpr_mod.stamp
	$(FC) $(FCFLAGS) -I. $< hexpr_mod.o $(LIBS) -o $@

test_completeness: test_completeness.f90 hexpr_mod.stamp
	$(FC) $(FCFLAGS) -I. $< hexpr_mod.o $(LIBS) -o $@

# ---- F77 tests (fixed form, use underscore aliases in libhexpr) --
# Compiled without -Wall/-Wextra to suppress style warnings on F77 code.
test_f77_smoke: test_f77_smoke.f
	$(FC) -O2 $< $(LIBS) -o $@

test_f77_completeness: test_f77_completeness.f
	$(FC) -O2 $< $(LIBS) -o $@

test_f77_eval: test_f77_eval.f
	$(FC) -O2 $< $(LIBS) -o $@

# ---- Run all tests, report PASS/FAIL ----------------------------
check: all
	@pass=0; fail=0; \
	for t in $(TESTS); do \
	  LD_LIBRARY_PATH=../../c:$$LD_LIBRARY_PATH ./$$t \
	    && pass=$$((pass+1)) \
	    || fail=$$((fail+1)); \
	done; \
	echo ""; \
	echo "Results: $$pass PASS, $$fail FAIL"; \
	[ $$fail -eq 0 ]

# ---- Clean ------------------------------------------------------
clean:
	rm -f *.o *.mod hexpr_mod.stamp $(TESTS)
