# =====================================================================
#  Makefile -- builds libwigner.a and libwigner.$(SOEXT).
# =====================================================================

CC      ?= gcc
# Warning flags. Default is warn-only so the build never fails on a future
# compiler's new warnings. To gate locally: make WARN='-Wall -Wextra -pedantic -Werror'
WARN    ?= -Wall -Wextra -pedantic
CFLAGS   = -O2 -std=c99 $(WARN)
CPPFLAGS = -Iinclude -Isrc
LDLIBS   = -lm

# --- platform shared-library conventions (Linux .so / macOS .dylib) ---
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SOEXT        := dylib
  SHARED       := -dynamiclib
  INSTALL_NAME  = -Wl,-install_name,@rpath/$(@F)
else
  SOEXT        := so
  SHARED       := -shared
  INSTALL_NAME  =
endif

LIB_SRCS = wigner_init.c wigner_3j.c wigner_6j.c wigner_9j.c
LIB_OBJS = $(LIB_SRCS:.c=.o)

LIB_STATIC = ../libwigner.a
LIB_SHARED = ../libwigner.$(SOEXT)

.PHONY: all clean

all: $(LIB_STATIC) $(LIB_SHARED)

$(LIB_STATIC): $(LIB_OBJS)
	$(AR) rcs $@ $^
	cp wigner.h ..

$(LIB_SHARED): $(LIB_OBJS)
	$(CC) $(SHARED) $(INSTALL_NAME) -o $@ $^ $(LDLIBS)
	cp wigner.h ..

$(LIB_OBJS): wigner_internal.h wigner.h

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -fPIC -c -o $@ $<

clean:
	rm -f $(LIB_OBJS) $(LIB_STATIC) $(LIB_SHARED) ../wigner.h

