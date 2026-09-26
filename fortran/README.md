# hexpr -- Fortran binding

## Overview

This directory provides Fortran bindings for the **hexpr** C library, with
two usage paths that share the same `libhexpr.so`:

> macOS: the shared library is `libhexpr.dylib`, and `DYLD_LIBRARY_PATH`
> replaces `LD_LIBRARY_PATH`. The test Makefile bakes an rpath so the tests
> run with no environment variable set.

- **F90 path** (`iso_c_binding` module): `fortran/src/hexpr_mod.f90` exposes
  the full public API as standard Fortran 95 procedures and constants. Use
  `use hexpr_mod` in modern Fortran code. Handles are `type(c_ptr)`; the
  module manages C pointer lifetimes and converts C strings to fixed-length
  Fortran characters.

- **F77 path** (underscore aliases): `c/src/hexpr_fortran.c` exports one
  underscore-suffixed wrapper per public function (e.g. `hexpr_init_`,
  `hexpr_drt_create_`). F77 programs call these by name without any module,
  passing all arguments by address. No separate Fortran source file is needed
  beyond the caller.

For the broader project context, see the [project root README](../README.md).

## Requirements

- One of the following Fortran compilers:
  - gfortran 9 or later (GNU)
  - LLVM Flang 20 or later
  - Intel ifx (oneAPI 2024+) -- requires `source setvars.sh` before both
    compilation and running (see [Caveats](#caveats))
- GNU Make
- `libhexpr.so` -- built from `c/` (see [c/README.md](../c/README.md))

## Installation

There is no separate package manager step. See
[Building from source](#building-from-source) -- the Makefile compiles the
F90 module and links the test programs in one shot.

## Quick start

```fortran
program quick_start
  use hexpr_mod
  use, intrinsic :: iso_c_binding
  implicit none
  integer(c_int) :: status

  status = hexpr_init()
  if (status /= HEXPR_OK) stop 1

  ! ... see workflows below ...

  call hexpr_shutdown()
end program quick_start
```

### 1. SOCI / FOCI workflow

Build the DRT for a FOCI or SOCI space and read the full expression in
one shot (this example uses the `s_full_4o4e` reference case: 4 valence
electrons, singlet SOCI):

```fortran
type(c_ptr) :: drt

drt = hexpr_drt_create( &
    kind    = HEXPR_DRT_SOCI, &
    nve     = 4_c_int,        &
    spin_x2 = 0_c_int,        &
    ncore   = 0_c_int,        &
    nval    = 4_c_int,        &
    next    = 0_c_int)

if (.not. c_associated(drt)) then
  write(*,*) 'hexpr_drt_create failed: ', trim(hexpr_last_error())
  stop 1
end if

write(*,'(a,i0,a,i0,a,i0)') 'norb=',  hexpr_drt_norb(drt),   &
                             '  nnodes=', hexpr_drt_nnodes(drt), &
                             '  ncsfs=',  hexpr_drt_ncsfs(drt)
! Expected: norb=4  nnodes=14  ncsfs=20

call hexpr_drt_destroy(drt)
```

Compile and run (from the repository root):

```sh
gfortran -I fortran/src -J fortran/src \
    fortran/src/hexpr_mod.f90 quick_start.f90 \
    -L c -l hexpr -o quick_start
LD_LIBRARY_PATH=c ./quick_start
```

### 2. CSF set workflow (per-pair / block / subspace evaluation)

Pass a chosen CSF set directly and evaluate individual pairs, rectangular
blocks, or square subspaces. Each step code is one of `0=E`, `1=U`,
`2=D`, `3=F` (see "CSF step codes" below). All CSFs must reach the same
top node.

```fortran
type(c_ptr)               :: drt
integer(c_int8_t), target :: csfs(9)
integer(c_int)            :: n, status
integer(c_int32_t), allocatable :: ij(:), pqrs(:)
real(c_double),     allocatable :: coef(:)

! Three CSFs over 3 orbitals: F F E, F U D, F E F.
! All reach the same top node (nve=4, spin=0).
csfs = [ 3_c_int8_t, 3_c_int8_t, 0_c_int8_t, &
         3_c_int8_t, 1_c_int8_t, 2_c_int8_t, &
         3_c_int8_t, 0_c_int8_t, 3_c_int8_t ]

drt = hexpr_drt_from_csfset(norb=3_c_int, csfs=csfs, n_csfs=3_c_int)
if (.not. c_associated(drt)) then
  write(*,*) 'hexpr_drt_from_csfset failed: ', trim(hexpr_last_error())
  stop 1
end if

write(*,'(a,i0,a,i0)') 'sub-DRT: nnodes=', hexpr_drt_nnodes(drt), &
                        '  ncsfs=',         hexpr_drt_ncsfs(drt)

! Evaluate the (bra=0, ket=1) pair of this DRT.
! First, get an upper bound on the term count to size the buffers.
status = hexpr_pair_max_terms(drt, HEXPR_EXPR_FULL, n)
allocate(ij(n), pqrs(n), coef(n))

! Then evaluate;
status = hexpr_pair_eval(drt, bra=0_c_int, ket=1_c_int, &
                         which=HEXPR_EXPR_FULL, out_count=n, &
                         out_ij=ij, out_pqrs=pqrs, out_coef=coef)

write(*,'(a,i0,a)') 'pair_eval returned ', n, ' terms'

deallocate(ij, pqrs, coef)
call hexpr_drt_destroy(drt)
```

`hexpr_block_eval` and `hexpr_subspace_eval` work analogously, taking
arrays of CSF row indices and returning extra `bra_pos` / `ket_pos`
arrays that index into the input arrays (see the API reference below).

## Index base

`ij` and `pqrs` are **0-based**, matching the C binary API. `ij` is a 0-based
*packed triangle index* of the `(bra, ket)` CSF pair (not a row index) — decode
it before use:

```fortran
i = int((sqrt(8.0d0 * ij + 1.0d0) - 1.0d0) / 2.0d0)   ! 0-based bra, i >= j
j = ij - i * (i + 1) / 2                               ! 0-based ket
```

then add 1 for Fortran's 1-based arrays (`H(i+1, j+1)`). `pqrs` is 0-based too.
The text dump (`reference/*.expr.txt`) adds 1 to `ij` only (to match the Ruby
reference corpus); `pqrs` stays 0-based.

## CSF step codes

CSFs are represented as sequences of step codes, one per orbital. The
step codes correspond to Shavitt's step operators `d` and his
`(a, b, c)` representation; this binding's `0/1/2/3` are essentially
the same as the `d` values used in Shavitt's Distinct Row Table
formulation, with minor differences in convention.

| Code | Symbol | Meaning |
|------|--------|---------|
| `0` | `E` | empty orbital (unoccupied) |
| `1` | `U` | singly occupied, spin-up coupling |
| `2` | `D` | singly occupied, spin-down coupling |
| `3` | `F` | doubly occupied, singlet coupling |

A CSF is well-formed when (a) every step code is in `0..3` and (b) the
cumulative `(a, b, c)` walked along the sequence stays in the physical
region (all three coordinates nonnegative at every orbital). The library
validates these conditions; CSFs reaching a common top node are required
when constructing a sub-DRT via `hexpr_drt_from_csfset`.

## API reference

The F90 module `hexpr_mod` exports the following public symbols. For full
argument lists and return types, see the source:
[`fortran/src/hexpr_mod.f90`](src/hexpr_mod.f90).

**Constants**

| Name | Value | Notes |
|------|-------|-------|
| `HEXPR_OK` | 0 | Success status |
| `HEXPR_ERR_INVALID_ARG` ... `HEXPR_ERR_INTERNAL` | 1...99 | Error codes |
| `HEXPR_DRT_SOCI` / `HEXPR_DRT_FOCI` / `HEXPR_DRT_CSFSET` | 0 / 1 / 2 | DRT kind |
| `HEXPR_EXPR_FULL` / `HEXPR_EXPR_ONE` | 0 / 1 | Expression selector |
| `HEXPR_STEP_E/U/D/F` | 0 to 3 | CSF step codes |
| `HEXPR_VMAJOR` / `HEXPR_VMINOR` / `HEXPR_VPATCH` | 1 / 0 / 0 | Version (see Caveats) |

**Procedures**

| Name | Kind | Description |
|------|------|-------------|
| `hexpr_init()` | function | Initialize library; returns status |
| `hexpr_shutdown()` | subroutine | Release library resources |
| `hexpr_version_string()` | function | Version as `character(len=64)` |
| `hexpr_version_major/minor/patch()` | function | Version integers |
| `hexpr_last_error()` | function | Last error as `character(len=256)` |
| `hexpr_drt_create(kind, nve, spin_x2, ncore, nval, next)` | function | Build DRT from FOCI/SOCI parameters |
| `hexpr_drt_from_csfset(norb, csfs, n_csfs)` | function | Build (sub-)DRT from an explicit CSF set |
| `hexpr_drt_destroy(drt)` | subroutine | Free DRT; sets handle to `c_null_ptr` |
| `hexpr_drt_norb/nnodes/ncsfs/kind/ncore/nval/next/nve_state/spin_x2(drt)` | function | DRT scalar queries |
| `hexpr_drt_node_info(drt, idx, orb, nve, is)` | function | Node details by index |
| `hexpr_csf_index(drt, step, out_index)` | function | CSF lookup: step codes -> 0-based index; returns status |
| `hexpr_csf_steps(drt, index, out_step)` | function | CSF lookup: index -> `norb` step codes into caller array; returns status |
| `hexpr_step_parse(str, norb, out_step)` | function | Parse `"eudf..."` string into step codes; returns status |
| `hexpr_csf_ndets(drt, index, ms_x2, out_ndets)` | function | Determinants CSF `index` expands into at `ms_x2`; exact, not a bound; returns status |
| `hexpr_csf_to_dets(drt, index, ms_x2, out_count, alpha, beta, coef)` | function | Determinant expansion -- see [CSF to determinants](#csf-to-determinants) |
| `hexpr_drt_csfs(drt, buf)` | function | CSF step-code matrix as 2-D `pointer` (CSFSET: input order) |
| `hexpr_drt_walk_csfs(drt, buf)` | function | Canonical GUGA-walk CSF matrix as 2-D `pointer` (CSFSET: the superset/canonical view) |
| `hexpr_drt_walk_ncsfs(drt)` | function | Canonical CSF count; CSFSET: may differ from `ncsfs` (superset M≥N or dedup) |
| `hexpr_drt_show(drt)` | function | Print DRT to stdout; returns status |
| `hexpr_expr_build(drt)` | function | Build expression; returns `type(c_ptr)` |
| `hexpr_expr_destroy(expr)` | subroutine | Free expression; sets handle to `c_null_ptr` |
| `hexpr_expr_nterms_full/one(expr)` | function | Term counts (`integer(c_int64_t)`) |
| `hexpr_expr_read(expr, which, ij, pqrs, coef)` | function | Bulk read with explicit selector |
| `hexpr_expr_read_full(expr, ij, pqrs, coef)` | function | Bulk read, full expression |
| `hexpr_expr_read_one(expr, ij, pqrs, coef)` | function | Bulk read, 1-electron subset |
| `hexpr_expr_show_full(expr)` | function | Print full expression to stdout; returns status |
| `hexpr_expr_show_one(expr)` | function | Print 1-electron expression to stdout; returns status |
| `hexpr_pair_count(drt, bra, ket, which, out_count)` | function | Term count for one CSF pair |
| `hexpr_pair_max_terms(drt, which, out_max)` | function | Upper bound on terms for any single pair |
| `hexpr_pair_eval(drt, bra, ket, which, out_count, ij, pqrs, coef)` | function | Pair evaluation |
| `hexpr_block_count(drt, bra_idx, n_bra, ket_idx, n_ket, which, out_count)` | function | Term count for a rectangular block |
| `hexpr_block_max_terms(drt, n_bra, n_ket, which, out_max)` | function | Upper bound on terms for a block |
| `hexpr_block_eval(drt, bra_idx, n_bra, ket_idx, n_ket, which, out_count, bra_pos, ket_pos, ij, pqrs, coef)` | function | Block evaluation |
| `hexpr_subspace_count(drt, idx, n, which, out_count)` | function | Term count for a square subspace |
| `hexpr_subspace_max_terms(drt, n, which, out_max)` | function | Upper bound on terms for a subspace |
| `hexpr_subspace_eval(drt, idx, n, which, out_count, bra_pos, ket_pos, ij, pqrs, coef)` | function | Subspace evaluation |

**F77 path**: underscore-suffixed names (`hexpr_init_`, `hexpr_drt_create_`,
`hexpr_pair_eval_`, etc.) are callable from Fortran 77 programs by linking
against `libhexpr.so` directly, without `use hexpr_mod`. All arguments are
passed by address. The scalar DRT getters (`hexpr_drt_norb_`,
`hexpr_drt_ncore_`, `hexpr_drt_nval_`, `hexpr_drt_next_`,
`hexpr_drt_nve_state_`, `hexpr_drt_spin_x2_`, and the rest) return their
integer result through a trailing `int *out` argument, consistent with all
other arguments being passed by address.
Source: [`c/src/hexpr_fortran.c`](../c/src/hexpr_fortran.c).

### Notes on the block / subspace return values

`hexpr_block_eval` and `hexpr_subspace_eval` write extra `bra_pos` and
`ket_pos` arrays alongside `(ij, pqrs, coef)`. For each output term,
`bra_pos(k)` and `ket_pos(k)` give the **position within the input index
arrays** (0-based), not the raw DRT CSF index. For example, if you call
`hexpr_block_eval` with `bra_idx = [5, 9, 12]` and `ket_idx = [3, 4]`,
then `bra_pos(k)` is in `{0, 1, 2}` and `ket_pos(k)` is in `{0, 1}`.
To recover the raw CSF index: `bra_idx(bra_pos(k) + 1)` (Fortran 1-based
indexing applied to the input array).

## Note: internal sub-DRT mechanism

`hexpr_drt_from_csfset` wraps the C function of the same name, which
constructs a DRT whose nodes and arcs are restricted to those reachable
from the supplied CSF set. This is the same internal mechanism used by
the C library's per-pair / block / subspace evaluation APIs
(`hexpr_pair_eval` and friends), exposed here for callers who want to
construct a sub-DRT directly for prototyping, custom analyses, or working
with externally-supplied CSF subsets.

The returned `type(c_ptr)` satisfies all standard DRT accessors (`norb`,
`nnodes`, `ncsfs`, etc.) and can be passed to any procedure that takes a
DRT handle, including `hexpr_expr_build` and the per-pair API. Its kind
is `HEXPR_DRT_CSFSET`.

**CSF ordering (sub-DRT `csfs`):** for a CSFSET sub-DRT from
`hexpr_drt_from_csfset`, `hexpr_drt_csfs` / `hexpr_drt_ncsfs` /
`hexpr_csf_index` / `hexpr_csf_steps` speak the **input set in input
order**: `hexpr_drt_csfs` returns exactly the CSFs you supplied, in the
order you supplied them, and `hexpr_drt_ncsfs` is that input count. The
`bra_pos`/`ket_pos` returned by `hexpr_block_eval`/`hexpr_subspace_eval`
index that same input array — consistent with `csfs`.

The internal DRT may expand to a canonical **superset** (more rows than the
input, when input CSFs share nodes) or dedup (fewer). That canonical view —
in canonical GUGA order (e<u<d<f), the form used internally for evaluation —
is available via **`hexpr_drt_walk_csfs` / `hexpr_drt_walk_ncsfs`**.

Validation of the input CSFs (step codes in 0..3, common top node) is
performed on the C side; failure returns `c_null_ptr` and sets
`hexpr_last_error()`.

## CSF to determinants

A CSF is a fixed linear combination of Slater determinants.
`hexpr_csf_ndets` reports how many determinants a CSF expands into at a
chosen `Ms`, and `hexpr_csf_to_dets` writes those determinants with their
coefficients into caller-supplied arrays. This is the seam to any
determinant-basis code: it lets a CSF-basis quantity be re-expressed in
the determinant basis, or a determinant-basis wavefunction be projected
onto CSFs.

```fortran
use hexpr_mod
use, intrinsic :: iso_c_binding
type(c_ptr)                             :: drt
integer(c_int)                          :: st, idx, norb, nwords, nd, p
integer(c_int8_t), target               :: parsed(4)
integer(c_int64_t), allocatable, target :: alpha(:,:), beta(:,:)
real(c_double),     allocatable, target :: coef(:)
logical                                 :: has_alpha

drt  = hexpr_drt_create(HEXPR_DRT_SOCI, 2_c_int, 0_c_int, &
                        1_c_int, 1_c_int, 2_c_int)
norb = hexpr_drt_norb(drt)
st   = hexpr_step_parse("udee", norb, parsed)
st   = hexpr_csf_index(drt, parsed, idx)

st = hexpr_csf_ndets(drt, idx, 0_c_int, nd)     ! ms_x2 = 0; here nd = 2
nwords = (norb + 63) / 64
allocate(alpha(nwords, nd), beta(nwords, nd), coef(nd))
st = hexpr_csf_to_dets(drt, idx, 0_c_int, nd, alpha, beta, coef)

p = 0                                            ! spatial orbital, 0-based
has_alpha = btest(alpha(ishft(p, -6) + 1, 1), iand(p, 63))
```

`ms_x2` is `2*Ms`. Any `|Ms| <= S` with the same parity as `spin_x2` is
accepted, negative values included -- `Ms` is not assumed to equal `S`.

`hexpr_csf_ndets` is **exact, not an upper bound**: it is
`C(n_open, n_alpha_open)`, so the arrays can be allocated from it with no
slack. This is unlike `hexpr_pair_max_terms` and friends, which are
bounds. Buffer lengths are a caller contract and are not diagnosed, as
with `hexpr_expr_read` and `hexpr_csf_steps`.

The `index` argument follows the same convention as `hexpr_csf_steps`: the
input position, in `[0, hexpr_drt_ncsfs)`, 0-based with no base conversion.
So do `hexpr_pair_eval` / `hexpr_block_eval` / `hexpr_subspace_eval` --
there is one index space across the whole API. For a CSFSET sub-DRT the
library maps input positions onto the canonical walk internally, so you
never supply a `hexpr_drt_walk_csfs` position.

### Determinant encoding -- and the one hazard Fortran has

`alpha(:, k)` and `beta(:, k)` are occupation bit masks for determinant
`k`, split into 64-bit words: bit `mod(p, 64)` of word `p / 64 + 1` is
spatial orbital `p`, with `nwords = (norb + 63) / 64`.

The C buffers are `out_alpha[k*nwords + w]` (row-major, word index
fastest). Fortran's *first* subscript is fastest, so the shape is
`(nwords, ndets)` and the same bytes read as `alpha(w, k)` -- no copy and
no transpose, the identical reinterpretation `hexpr_drt_csfs` already
uses.

**`iso_c_binding` has no `c_uint64_t`.** The C side's `uint64_t` therefore
arrives as *signed* `integer(c_int64_t)`, so a determinant occupying
spatial orbital 63 sets the sign bit and the mask reads as a negative
number. Test occupation with `BTEST`, which works on the bit pattern.
Never with a magnitude comparison, a decimal `write`, or a use as an array
subscript. This hazard exists in neither C nor the Python, Julia or Ruby
bindings, all of which have an unsigned or arbitrary-precision type.

### Conventions

These are fixed by the C API and derived in `c/src/csf/csf2det.c`. A
caller comparing against another program must match them, or the
coefficients will not agree.

- **Spin coupling**: genealogical, in the Eq. (5) order -- the new
  one-orbital state couples on the **left** of the accumulated partial
  CSF. A CSF expressed in another coupling scheme will not reproduce these
  coefficients.
- **Clebsch-Gordan phase**: Condon-Shortley.
- **Operator order**: all alpha in ascending spatial orbital, then all
  beta in ascending spatial orbital. An interleaved-by-orbital convention
  (used by some qubit mappings) differs from this by a
  determinant-dependent permutation sign.
- **Determinant order**: lexicographically ascending in the set of
  open-shell positions carrying alpha.
- **Normalization**: the coefficients of any one CSF square-sum to 1 at
  every `Ms`.
- **Zero coefficients** are emitted, not dropped, so the count stays exact
  and no threshold enters the API.

The **overall sign** of a CSF is a convention, not a measured quantity: it
is whatever the coupling order and the operator order produce. Relative
signs -- within one CSF, and between CSFs of one DRT -- are meaningful,
and are verified against PySCF by `validation/csf2det_vs_pyscf.py`. Both
`<I|J>` and `U H_det U^T` are invariant under flipping the sign of a whole
CSF, so the overall sign is not something to assert against.

### One shape, five languages

All five bindings expose the same `(nwords, ndets)` buffer and the same
conventions: the **data shape follows the C contract** rather than each
language's idiom. (Python, Julia and Ruby index it as `(ndets, nwords)` --
the same bytes, with their row-major subscript order.) Only error
reporting is adapted per language -- a status return here and in C,
exceptions in Python, `error()` in Julia, `raise` in Ruby.

**F77 path**: `hexpr_csf_ndets_` and `hexpr_csf_to_dets_` take the same
arguments by address, with a trailing `int *out_status`. The masks are
`integer*8` there, with the same signedness caveat.

## Note: no high-level "step-code in, expression out" wrappers

The Python, Julia, and Ruby bindings provide high-level functions
(`csf_pair_eval`, `csf_block_eval`, `csf_subspace_eval`) that take CSF
step-code sequences directly, build a transient sub-DRT internally,
evaluate, and discard the sub-DRT. Those languages have garbage
collection, so the transient DRT can be created and dropped invisibly
to the user.

Fortran does not have GC, so this style does not fit the language
convention. Fortran users call `hexpr_drt_from_csfset` explicitly, then
the per-pair / block / subspace API, then `hexpr_drt_destroy`. The
explicit DRT lifecycle matches Fortran's general memory-management
contract and lets the caller reuse the DRT across many evaluations.

This concerns only the `csf_*_eval` family. It is **not** a general
statement that Fortran gets less of the API: `hexpr_csf_ndets` and
`hexpr_csf_to_dets`, for instance, take an existing DRT and are present
here exactly as in the other four bindings.

## Building from source

The `fortran/tests/` Makefile compiles the F90 module and runs all test
programs. The default compiler is gfortran:

```sh
make -C fortran/tests check                           # gfortran (default)
make -C fortran/tests FC=flang clean check            # LLVM Flang
source /opt/intel/oneapi/setvars.sh \
  && make -C fortran/tests FC=ifx clean check         # Intel ifx
```

The Makefile produces the compiled module file (`hexpr_mod.mod`) in the
`fortran/tests/` working directory.

## Running tests

The test command above builds and runs all tests in one step. Test
programs cover:

| Program | Coverage |
|---------|----------|
| `test_smoke` | `hexpr_init`, version query, `hexpr_shutdown` |
| `test_drt` | DRT construction and scalar accessors for `s_full_4o4e` |
| `test_expr` | Expression build, term counts, `read_full`/`read_one` |
| `test_csfset` | `hexpr_drt_from_csfset`: round-trip and 2-CSF sub-DRT |
| `test_pair` | Per-pair / block / subspace API including aggregate vs full expression |
| `test_completeness` | F90 binding-completeness wrappers (`hexpr_csf_index`, `hexpr_csf_steps`, `hexpr_step_parse`, the `*_show*` procedures) and the determinant expansion |
| `test_f77_smoke` | F77 path: `hexpr_init_` / `hexpr_shutdown_` via aliases |
| `test_f77_completeness` | F77 binding-completeness aliases and the determinant expansion |
| `test_f77_eval` | F77 path: per-pair evaluation against a reference CSF set |

The determinant-expansion checks in `test_completeness` and
`test_f77_completeness` assert determinant counts, bit masks, magnitudes
and *relative* signs. They never assert the sign of an individual
coefficient: the overall sign of a CSF is a convention (see
[CSF to determinants](#csf-to-determinants)). `test_completeness` reads
the occupation back with the `BTEST` expression the module documents,
which is the form that still works at `norb > 64`.

All tests pass on gfortran, flang, and ifx.

## See also

- Project root: [`../README.md`](../README.md)
- C library: [`../c/README.md`](../c/README.md)
- Public API header: [`../c/include/hexpr.h`](../c/include/hexpr.h)
- F90 module source: [`fortran/src/hexpr_mod.f90`](src/hexpr_mod.f90)
- Index encoding: [`../docs/PQRS_INDEX_SPEC.md`](../docs/PQRS_INDEX_SPEC.md)
- Python binding: [`../python/README.md`](../python/README.md)
- Julia binding: [`../julia/README.md`](../julia/README.md)
- Ruby binding: [`../ruby/README.md`](../ruby/README.md)
- Methodology paper: [DOI placeholder until preprint is up]

## Caveats

**`HEXPR_VMAJOR` / `HEXPR_VMINOR` / `HEXPR_VPATCH` (not `HEXPR_VERSION_MAJOR`)**
Fortran is case-insensitive. The name `HEXPR_VERSION_MAJOR` would clash with
the runtime function `hexpr_version_major()`. The module therefore uses the
short-form constants `HEXPR_VMAJOR`, `HEXPR_VMINOR`, `HEXPR_VPATCH`.

**`c_loc` requires contiguous arrays (F2003)**
`hexpr_expr_read`, `hexpr_drt_from_csfset`, and the per-pair / block /
subspace evaluators use `c_loc` to pass assumed-shape array addresses to
C. This requires the actual arguments to be contiguous in memory.
Allocatable arrays always satisfy this; Fortran pointer arrays may need
the `CONTIGUOUS` attribute. There is no static check in F95-compatible
mode.

**`INTENT(OUT)` not allowed on assumed-size array dummies in `bind(C)` interfaces**
A Fortran 95 constraint: assumed-size array dummies in `bind(C)`
interface blocks cannot carry `INTENT(OUT)`. The bulk-read interface in the
module passes the output arrays as `type(c_ptr), value` to work around this.

**flang prints "IEEE arithmetic exceptions signaled: INEXACT"**
After each `STOP` statement in a test program, flang 20 prints this warning
to stderr. It is a harmless artifact of floating-point rounding in standard
library finalization; it cannot be suppressed without changing visible program
behavior.

**gfortran prints "STOP 0" on successful exit**
gfortran's default behavior is to print `STOP 0` to stderr after each
successful `STOP 0` statement. This is expected and harmless.

**ifx requires `source setvars.sh`**
Intel ifx needs `libimf.so` and related runtime libraries from the oneAPI
installation. These are not on the default `LD_LIBRARY_PATH`. Run
`source /opt/intel/oneapi/setvars.sh` in your shell before both compiling
and running ifx-built binaries.
