# Wigner 3j / 6j / 9j coefficient library

Original Fortran implementation by F. Sasaki, refactored into modern C
while preserving bit-exact numerical agreement with the original Fortran
output.

Copyright belongs to F. Sasaki, separate from the main project copyright
(H. Honda, T. Noro), under the same MIT license terms -- see
[LICENSE](LICENSE).

## Public C API (with bonus Ruby extension)

`wigner.h` declares `wigner_init`, `wigner_delta`, `wigner_f3l`,
`wigner_6j`, and `wigner_9j` in `extern "C"`-safe form. The Ruby extension
files (`wigner_ruby.h`, `wigner_ruby.c`, `extconf.rb`) are co-located in
this directory as a bonus, with a separate build flow (see "Ruby
extension" below).

## Files
wigner.h            public header (the only file callers need to include)
wigner_internal.h   internal header (table structs, wig_phase, etc.)
wigner_init.c       table initialisation
wigner_3j.c         wigner_delta, wigner_f3l
wigner_6j.c         wigner_6j
wigner_9j.c         wigner_9j
wigner_ruby.h       Ruby extension header (bonus)
wigner_ruby.c       Ruby extension implementation (bonus)
extconf.rb          Ruby extension build config (bonus)
test_wigner.c       exhaustive comparison test against the original Fortran
                    (f2c) code -- see "Exhaustive test" below; the sources it
                    compares against are NOT in this repository
example.c           example of calling the library from C
Makefile            builds libwigner.a / libwigner.so
LICENSE             F. Sasaki copyright (separate from main project)

## Public API

```c
#include "wigner.h"

void   wigner_init(void);

double wigner_delta(int two_ja_plus_1, int two_jb_plus_1, int two_jc_plus_1);
double wigner_f3l  (int two_ja_plus_1, int two_jb_plus_1, int two_jc_plus_1);

double wigner_6j(int two_j1_plus_1, int two_j2_plus_1, int two_j3_plus_1,
                 int two_j4_plus_1, int two_j5_plus_1, int two_j6_plus_1);

double wigner_9j(int two_ja_plus_1, int two_jb_plus_1, int two_jc_plus_1,
                 int two_jd_plus_1, int two_je_plus_1, int two_jv_plus_1,
                 int two_jg_plus_1, int two_jh_plus_1, int two_ji_plus_1);
```

Arguments follow the same convention as the original code: each is the integer
value of `(2j+1)`. Half-integer `j` is handled via odd/even parity.

Call `wigner_init()` once before any other routine.

## Building

This library is built automatically as part of the main hexpr build
(`make -C c`). To build standalone:

```sh
make        # produces libwigner.a and libwigner.so
make clean
```

The Makefile builds only the C library. `test_wigner.c` has no rule here;
see "Exhaustive test" below. The Ruby extension is not part of this
Makefile and is built separately (see below).

To link from a C program:

```sh
cc -I/path/to/c/wig  myprog.c  /path/to/c/libwigner.a  -lm
# or
cc -I/path/to/c/wig  myprog.c  -L/path/to/c  -lwigner  -lm
```

## Ruby extension (bonus)

The Ruby extension is bundled as a bonus and built via Ruby's standard
extension workflow (separate from the main `Makefile`):

```sh
ruby extconf.rb
make
# produces Wigcoef.so (or .bundle)
```

Module and function names are compatible with the original Ruby wrapper
(`Wigcoef.init`, `.delta`, `.f3l`, `.wig6j`, `.wig9x`).
The 9j symbol is also accessible as `.wig9j`; `.wig9x` is kept as an alias.

## Note for source modifiers

One subtle bug to be aware of when modifying the table initialisation
code. The original `setwig_` contains:

```c
fasq[(k << 1) + 38] = den[k - 1];   // for k = 1..57 (1-based)
```

When converting to 0-based indexing the **offset must be 40, not 38** because
at `k = 1` the target index is `2*1 + 38 = 40`. Writing 38 instead leaves
`delta`'s index arithmetic unchanged but shifts the table values, so forbidden
triangles that should return 0 return non-zero instead. Because `wigner_6j`
uses `delta == 0` to detect forbidden triangles, this silently corrupts 6j and
9j values for combinations that look geometrically valid but are actually
forbidden.

The exhaustive test in `test_wigner.c` catches this bug, and anyone
modifying this code should run it.

## Exhaustive test

`test_wigner.c` sweeps every entry point over a range of `(2j+1)` values
and compares each result against the original f2c-translated routines
(`setwig_`, `delta_`, `f3l_`, `xwig6j_`, `xwig9x_`), which it declares
`extern`. It is the only check that covers the failure mode described
above: the three spot values in `c/tests/test_wigner.c` do not.

**Those f2c sources are not in this repository, so this file cannot be
compiled as it stands.** The `Makefile` has no `test` target. To run it,
supply the original translated sources and link them alongside this file
and the library:

```sh
cc -I. -o test_wigner test_wigner.c <f2c objects> ../libwigner.a -lm
```

Replacing the f2c side with a checked-in table of expected values would
make this runnable from a clean checkout. That has not been done.

## Limitations

- Only `(2j+1) <= 56` is supported. Exceeding this causes `wigner_6j` to print
  a warning to stderr and return 0 (preserving the original behaviour).
- The table filled by `wigner_init()` is in static storage. Reading from
  multiple threads is safe; `wigner_init()` itself is not re-entrant. Call it
  once at program startup.
- `wigner_delta` returns only the product of three `sqrt((+/-)!)` terms -- the
  `(a+b+c+1)!` denominator is applied on the `wigner_f3l` side. This differs
  from the standard triangle coefficient definition; keep this in mind when
  calling `wigner_delta` directly.
