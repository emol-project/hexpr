# Default compiler is gcc. Override with `make CC=clang` or
# `make CC=icx`. For icx, source the Intel oneAPI environment
# first (`source /opt/intel/oneapi/setvars.sh`); this sets
# LD_LIBRARY_PATH for the Intel compiler runtime libraries.
CC      ?= gcc
# Warning flags. Default is warn-only so the build never fails on a future
# compiler's new warnings. To gate locally: make WARN='-Wall -Wextra -pedantic -Werror'
WARN    ?= -Wall -Wextra -pedantic

# -fvisibility=hidden: nothing is exported from libhexpr.so unless a header
# says so. include/hexpr.h and src/hexpr_fortran.h mark their declarations
# with a visibility pragma; everything else becomes internal to the library.
CFLAGS   = $(WARN) -O2 -g -fPIC -std=c11 -fvisibility=hidden
CPPFLAGS = -Iinclude -I. -Iwig -Isrc/common -Isrc/loopgen -Isrc/csf
LDLIBS   = -lm

# --- platform shared-library conventions (Linux .so / macOS .dylib) ---
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SOEXT          := dylib
  SHARED         := -dynamiclib
  INSTALL_NAME    = -Wl,-install_name,@rpath/$(@F)
  WHOLE_ARCHIVE   = -Wl,-force_load,$(WIG_ARCHIVE)
  RPATH_BIN      := -Wl,-rpath,@loader_path
  # No Darwin equivalent of --exclude-libs. The wigner symbols therefore
  # stay visible in libhexpr.dylib; the c/src ones are hidden as on Linux.
  EXCLUDE_WIG    :=
else
  SOEXT          := so
  SHARED         := -shared
  INSTALL_NAME    =
  WHOLE_ARCHIVE   = -Wl,--whole-archive $(WIG_ARCHIVE) -Wl,--no-whole-archive
  RPATH_BIN      :=
  EXCLUDE_WIG    := -Wl,--exclude-libs,libwigner.a
endif

COMMON_SRCS = src/common/fmt.c src/common/darr.c src/common/dhash.c src/common/sort.c
COMMON_OBJS = $(COMMON_SRCS:.c=.o)

LOOPGEN_SRCS = src/loopgen/drt.c src/loopgen/expr_tree.c src/loopgen/expr.c
LOOPGEN_OBJS = $(LOOPGEN_SRCS:.c=.o)

CSF_SRCS = src/csf/csf_lookup.c src/csf/csf_pair_core.c src/csf/csf_pair.c src/csf/csf2det.c
CSF_OBJS = $(CSF_SRCS:.c=.o)

# Core API sources.
API_SRCS_CORE = src/hexpr_api.c src/hexpr_fortran.c
API_OBJS_CORE = $(API_SRCS_CORE:.c=.o)

WIG_SRCS = wig/wigner_init.c wig/wigner_3j.c wig/wigner_6j.c wig/wigner_9j.c
WIG_OBJS = $(WIG_SRCS:.c=.o)
WIG_ARCHIVE = libwigner.a

LIBHEXPR      = libhexpr.$(SOEXT)
LIBHEXPR_A    = libhexpr.a
RUN_REF       = run_reference

.PHONY: all clean

# Default: libhexpr.so + run_reference.
all: $(LIBHEXPR) $(LIBHEXPR_A) $(RUN_REF)

$(WIG_ARCHIVE): $(WIG_SRCS) wig/wigner.h wig/wigner_internal.h
	$(MAKE) -C wig

src/common/%.o: src/common/%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

src/loopgen/%.o: src/loopgen/%.c src/loopgen/drt.h src/loopgen/expr_internal.h src/loopgen/expr.h src/loopgen/brooks_table.h  src/loopgen/lk_table.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

src/csf/%.o: src/csf/%.c src/csf/csf_lookup.h src/csf/csf_pair_core.h src/loopgen/expr_internal.h src/loopgen/brooks_table.h  src/loopgen/lk_table.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

src/hexpr_api.o: src/hexpr_api.c include/hexpr.h wig/wigner.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

src/hexpr_fortran.o: src/hexpr_fortran.c include/hexpr.h src/hexpr_fortran.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(LIBHEXPR): $(COMMON_OBJS) $(LOOPGEN_OBJS) $(CSF_OBJS) $(API_OBJS_CORE) $(WIG_ARCHIVE)
	$(CC) $(SHARED) -fPIC $(INSTALL_NAME) -o $@ \
	    $(COMMON_OBJS) $(LOOPGEN_OBJS) $(CSF_OBJS) $(API_OBJS_CORE) \
	    $(EXCLUDE_WIG) $(WHOLE_ARCHIVE) $(LDLIBS)

# Static archive over the same objects as the shared library, with the
# wigner objects folded in so it is self-contained (-lm is still needed by
# the caller). It exists because the shared library hides every function
# that include/hexpr.h does not declare, while some programs in c/tests
# deliberately exercise internal entry points; those link this instead.
# Depending on $(WIG_ARCHIVE) is what guarantees $(WIG_OBJS) have been
# built by the wig sub-make.
$(LIBHEXPR_A): $(COMMON_OBJS) $(LOOPGEN_OBJS) $(CSF_OBJS) $(API_OBJS_CORE) $(WIG_ARCHIVE)
	rm -f $@
	$(AR) rcs $@ $(COMMON_OBJS) $(LOOPGEN_OBJS) $(CSF_OBJS) \
	    $(API_OBJS_CORE) $(WIG_OBJS)
# run_reference calls hexpr_drt_mk_csf_list and hexpr_print_csf, which the
# shared library does not export. It links the static archive, like the
# c/tests programs that reach past include/hexpr.h.
$(RUN_REF): src/run_reference.c $(LIBHEXPR_A)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LIBHEXPR_A) $(LDLIBS)

clean:
	rm -f $(COMMON_OBJS) $(LOOPGEN_OBJS) $(CSF_OBJS) \
	      $(API_OBJS_CORE) \
	      $(LIBHEXPR) $(LIBHEXPR_A) $(WIG_ARCHIVE) $(RUN_REF)
	$(MAKE) -C wig clean
