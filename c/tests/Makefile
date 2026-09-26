CC      ?= gcc
# Warning flags. Default is warn-only so the build never fails on a future
# compiler's new warnings. To gate locally: make WARN='-Wall -Wextra -pedantic -Werror'
WARN    ?= -Wall -Wextra -pedantic
CFLAGS   = $(WARN) -O2 -g -std=c11
CPPFLAGS = -Iutil -I../include -I../wig -I../src/common -I../src -I../src/csf

# --- platform conventions (rpath so tests run without DYLD_*/LD_LIBRARY_PATH) ---
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  RPATH_TEST := -Wl,-rpath,@loader_path/..
  DL_LIB     :=
  RT_LIB     :=
else
  RPATH_TEST :=
  DL_LIB     := -ldl
  RT_LIB     := -lrt
endif

LDFLAGS  = -L.. $(RPATH_TEST)
LDLIBS   = -lhexpr -lm

# libhexpr.so exports only what ../include/hexpr.h declares plus the F77
# wrappers. A few programs here deliberately call internal entry points
# (drt.h, fmt.h, expr_internal.h); they link the static archive, which is
# built from the same objects as the shared library and already contains
# the wigner objects. Named by path, because -lhexpr would find the shared
# library first.
LDLIBS_STATIC = ../libhexpr.a -lm

# test_wigner exercises libwigner directly and does not touch libhexpr.
LDLIBS_WIGNER = ../libwigner.a -lm

LIBDIR   = $(abspath ..)

UTIL_OBJS = util/compare.o

TESTS_CORE = test_compare test_wigner test_drt test_expr test_full \
             test_api test_fortran_symbols test_csf_lookup \
             test_pair_one_4o4e test_pair_full_4o4e test_block_subspace_4o4e \
             test_pair_aggregate test_max_terms test_drt_from_csfset \
             test_csfset_pair test_ncsfs_consistency test_lk_table test_lk_twophase \
             test_csf2det

.PHONY: all check check-slow clean

all: $(TESTS_CORE)

util/compare.o: util/compare.c util/compare.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

test_compare: test_compare.c $(UTIL_OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_wigner: test_wigner.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS_WIGNER)

test_drt: test_drt.c $(UTIL_OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS_STATIC)

test_expr: test_expr.c $(UTIL_OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_full: test_full.c $(UTIL_OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_api: test_api.c $(UTIL_OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) -lm

test_fortran_symbols: test_fortran_symbols.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(DL_LIB)

test_csf_lookup: test_csf_lookup.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_pair_one_4o4e: test_pair_one_4o4e.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_pair_full_4o4e: test_pair_full_4o4e.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_block_subspace_4o4e: test_block_subspace_4o4e.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_pair_aggregate: test_pair_aggregate.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_max_terms: test_max_terms.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_drt_from_csfset: test_drt_from_csfset.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS_STATIC)

test_csfset_pair: test_csfset_pair.c csfset_oracle_util.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ test_csfset_pair.c $(LDLIBS)

test_ncsfs_consistency: test_ncsfs_consistency.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

test_csf2det: test_csf2det.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Self-contained: exercises the POD lk_table.h directly, no libhexpr link.
test_lk_table: test_lk_table.c ../src/loopgen/lk_table.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -I../src/loopgen -o $@ test_lk_table.c

test_lk_twophase: test_lk_twophase.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS_STATIC)

# Regenerate the committed reference_csfset/ oracle (transparency only;
# the committed files ARE the oracle and the test does not need this).
gen_csfset_oracle: gen_csfset_oracle.c csfset_oracle_util.h
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ gen_csfset_oracle.c $(LDLIBS)

.PHONY: gen-csfset-oracle
gen-csfset-oracle: gen_csfset_oracle
	LD_LIBRARY_PATH=$(LIBDIR) ./gen_csfset_oracle

test_expr_stress: test_expr_stress.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(RT_LIB)

check: $(TESTS_CORE)
	@failures=0; \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_compare          || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_wigner           || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_drt              || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_expr             || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_full             || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_api              || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_fortran_symbols  || failures=$$((failures+1)); \
	sh ./check_fortran_coverage.sh                    || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_csf_lookup       || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_pair_one_4o4e    || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_pair_full_4o4e          || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_block_subspace_4o4e     || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_pair_aggregate s_full_4o4e t_full_4o4e || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_max_terms                              || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_drt_from_csfset                       || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_csfset_pair                           || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_ncsfs_consistency                     || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_csf2det                               || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_lk_table                              || failures=$$((failures+1)); \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_lk_twophase                           || failures=$$((failures+1)); \
	if [ $$failures -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "$$failures TEST(S) FAILED"; fi; \
	exit $$failures

.PHONY: check-slow
check-slow: $(TESTS_CORE)
	@failures=0; \
	LD_LIBRARY_PATH=$(LIBDIR) ./test_pair_aggregate \
	    s_full_4o4e t_full_4o4e s_full_6o6e t_full_6o6e || failures=$$((failures+1)); \
	if [ $$failures -eq 0 ]; then echo "ALL SLOW TESTS PASSED"; else echo "$$failures TEST(S) FAILED"; fi; \
	exit $$failures

.PHONY: stress
stress: test_expr_stress
	LD_LIBRARY_PATH=$(LIBDIR) ./test_expr_stress

clean:
	rm -f $(UTIL_OBJS) $(TESTS_CORE)
	rm -rf tmp/
