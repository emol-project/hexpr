# hexpr -- C library

## Overview

hexpr is a C library that generates expressions of the CSF-based
Hamiltonian matrix for quantum chemistry CI calculations. Input CSFs
can be supplied directly, or the full CSF set of a first-order CI
(FOCI) or second-order CI (SOCI) space can be specified by its orbital
classification and hexpr will generate the corresponding CSFs internally.

The matrix element formulation follows F. Sasaki's tensor recoupling and reduction technique,
originally developed for atomic systems in LS coupling
and here applied to the molecular case without assuming spatial symmetry.
Internally, hexpr uses the Distinct Row Table (DRT) representation and
Brooks-style loop traversal.

The library accepts either an explicit list of spin-coupled CSFs, or an
orbital classification (number of core, valence, and external orbitals;
number of valence electrons; spin multiplicity), and returns the complete
set of Hamiltonian matrix elements in compressed form: parallel
integer arrays `(ij, pqrs, coef)` encoding CSF-pair indices, orbital indices,
and coefficients. Consumers combine these with one- and two-electron MO
integrals to assemble and diagonalize the CI Hamiltonian.

Language bindings (Python, Fortran, Julia, Ruby) are provided in sibling
directories. See the [project root README](../README.md) for an overview.

## Requirements

- C11 compiler (GCC 9+ or Clang 11+ recommended)
- GNU Make
- libm (standard; always available)

## Building

> macOS: the build produces `libhexpr.dylib` (not `.so`), and runtime uses
> `DYLD_LIBRARY_PATH` in place of `LD_LIBRARY_PATH`. The Makefile selects the
> right form automatically.

```sh
# from the repository root -- builds libhexpr.so and the run_reference driver
make -C c
```

Produces:

- `c/libhexpr.so` -- the shared library
- `c/run_reference` -- standalone CLI driver
- `c/libhexpr.a` -- static archive of the same objects, linked by
  `run_reference` and by the C test suites that call internal functions

To rebuild cleanly:

```sh
make -C c clean && make -C c
```

## API reference

The stable public API is declared in [`c/include/hexpr.h`](include/hexpr.h).
That header is the source of truth. Internal headers under `src/loopgen/`,
`src/csf/`, `src/common/`, and `wig/` are not part of the public contract.

This is enforced at link time. The library is compiled with
`-fvisibility=hidden`, so `libhexpr.so` exports only the functions
declared in `hexpr.h`, their F77 wrappers, and `hexpr_fortran_stdout`
(used by the Fortran module). A function declared only in an internal
header cannot be linked against the shared library.
`c/tests/check_fortran_coverage.sh`, run by `make -C c/tests check`,
compares the exported symbol table with both headers. On macOS the
bundled Wigner routines remain exported; see `CHANGELOG.md`.

Callers should:

1. `#include "hexpr.h"`
2. Call `hexpr_init()` once before any other library call
3. Link against `libhexpr.so` (`-L<path> -lhexpr -lm`)

Key types and functions at a glance:

| Symbol | Description |
|--------|-------------|
| `hexpr_init()` / `hexpr_shutdown()` | Library lifecycle |
| `hexpr_drt_create(kind, nve, spin_x2, ncore, nval, next)` | Build a DRT |
| `hexpr_drt_destroy(drt)` | Free a DRT |
| `hexpr_drt_norb/nnodes/ncsfs(drt)` | DRT scalar queries |
| `hexpr_drt_csfs(drt)` | CSF list as packed step-code bytes (CSFSET: input order) |
| `hexpr_drt_walk_csfs(drt)` / `hexpr_drt_walk_ncsfs(drt)` | Canonical GUGA-walk CSF list + count (CSFSET: the superset/canonical view) |
| `hexpr_expr_build(drt)` | Build the energy expression |
| `hexpr_expr_destroy(expr)` | Free an expression |
| `hexpr_expr_nterms_full/one(expr)` | Term counts |
| `hexpr_expr_read(expr, which, ij, pqrs, coef)` | Bulk read into caller buffers |
| `hexpr_pair_count/eval(drt, bra, ket, ...)` | Per-pair expression evaluation |
| `hexpr_block_count/eval(drt, bras, kets, ...)` | Block evaluation (rectangular CSF set) |
| `hexpr_subspace_count/eval(drt, indices, ...)` | Subspace evaluation (square CSF set) |
| `hexpr_pair/block/subspace_max_terms(...)` | Upper-bound for output buffer sizing |
| `hexpr_csf_index(drt, step, &index)` | CSF lookup: step codes -> index |
| `hexpr_csf_steps(drt, index, step)` | CSF lookup: index -> step codes |
| `hexpr_step_parse(str, norb, step)` | Parse "eudf..." string to bytes |
| `hexpr_csf_ndets(drt, index, ms_x2, &ndets)` | Determinant count for one CSF at a given Ms; exact, not a bound |
| `hexpr_csf_to_dets(drt, index, ms_x2, &count, alpha, beta, coef)` | CSF -> Slater determinant expansion |
| `hexpr_last_error()` | Thread-local error string |

## Examples

### CLI driver

`c/run_reference` runs all 14 built-in reference cases and writes three text
files per case (DRT table, full expression, 1-electron expression):

```sh
# All 14 cases -> ./reference_c/
LD_LIBRARY_PATH=c c/run_reference

# Selected cases, custom output directory
LD_LIBRARY_PATH=c c/run_reference -o /tmp/out s_full_4o4e t_full_4o4e
```

Per-case progress (case header, file sizes, elapsed time) goes to stderr;
stdout is silent so output can be redirected freely.

Source: [`c/src/run_reference.c`](src/run_reference.c)

### Wigner library example

The Wigner 3j/6j/9j sub-library in `c/wig/` has its own standalone example:

```sh
cd c/wig
make
./example
```

Source: [`c/wig/example.c`](wig/example.c)

## Library structure
c/
+-- include/hexpr.h          Public API header (v1.0.0, frozen)
+-- src/
|   +-- common/              Shared utilities: fmt, darr, dhash, sort
|   +-- loopgen/             DRT construction + Brooks expression
|   |                        (drt.c, expr.c, expr_tree.c)
|   +-- csf/                 CSF-pair API and CSF -> determinant expansion
|   |                        (csf_pair.c, csf_pair_core.c, csf_lookup.c,
|   |                         csf2det.c)
|   +-- hexpr_api.c          Version, hexpr_init/shutdown, hexpr_last_error
|   +-- hexpr_fortran.c      Underscore-suffixed F77 wrappers
|   `-- run_reference.c      CLI driver
`-- wig/                     Wigner 3j/6j/9j sub-library (F. Sasaki)

## Testing

```sh
make -C c             # build library and driver first
make -C c/tests check # run the core test suite
```

`make check` runs the core test suite, in the order defined by
`TESTS_CORE` in `c/tests/Makefile`.

Expected output ends with `ALL TESTS PASSED`.

Core test descriptions:

| Suite | What it covers |
|-------|----------------|
| `test_compare` | File-comparison utility self-test |
| `test_wigner` | Wigner 3j/6j/9j coefficient correctness |
| `test_drt` | DRT construction for 12 of the 14 reference cases (all but `*_mid`) |
| `test_expr` | Expression generation for 12 of the 14 reference cases (all but `*_mid`) |
| `test_full` | End-to-end driver vs. reference files (14 cases, 42 files) |
| `test_api` | Every public function in `hexpr.h` |
| `test_fortran_symbols` | Every F77 underscore-suffixed symbol declared in `src/hexpr_fortran.h` is present in libhexpr.so |
| `test_csf_lookup` | `hexpr_csf_index`, `hexpr_csf_steps`, `hexpr_step_parse` |
| `test_pair_one_4o4e` | Per-pair 1e-only expression evaluation |
| `test_pair_full_4o4e` | Per-pair full expression evaluation |
| `test_block_subspace_4o4e` | `block_eval` / `subspace_eval` |
| `test_pair_aggregate` | Pair-aggregate equals `expr_build` |
| `test_max_terms` | `hexpr_pair/block/subspace_max_terms` upper bounds |
| `test_drt_from_csfset` | DRT construction from explicit CSF set |
| `test_csfset_pair` | CSFSET evaluation against the frozen `reference_csfset/` oracle, including the input-order contract |
| `test_ncsfs_consistency` | `ncsfs` vs `walk_ncsfs`, and rejection of target states with no valid walk |
| `test_lk_table` | `L(k)` lookup table: mixed-radix packing, dedup, overfill and integrity aborts |
| `test_lk_twophase` | Two-phase `L(k)`: whole-path search must be lookup-only (no misses, no recomputation) |
| `test_csf2det` | `hexpr_csf_ndets` / `hexpr_csf_to_dets`: determinant counts, orthonormality, Sz and S^2 across every Ms, CG closed forms, relative signs, and multi-word (`norb > 64`) determinants |

Three further checks on the determinant expansion live outside this
suite. `validation/csf2det_vs_pyscf.py` compares the whole
determinant-basis matrix against PySCF for 14 cases -- the only check
that constrains the *sign* convention against anything outside hexpr,
since counts, magnitudes, norms, `<I|J>` and `S^2` are all invariant
under flipping the sign of a CSF. The Python, Julia and Ruby bindings
each carry their own determinant-expansion tests, and the Fortran ones
are in `fortran/tests/`.

### Reference data for full regression coverage

Four of the largest reference files are not bundled with the
repository to keep clone size manageable:

- `reference/s_foci_med.expr.txt`
- `reference/s_soci_med.expr.txt`
- `reference/t_foci_med.expr.txt`
- `reference/t_soci_med.expr.txt`

Without these files, the test suite still reports `ALL TESTS
PASSED`; the affected comparisons print `SKIP` instead of `PASS`.
That is enough to verify a build is working.

To run the full regression with all four `*_med.expr.txt` files
in place, download the complete reference set:

```sh
wget https://github.com/emol-project/hexpr/releases/download/v1.0.0/reference.tar.gz
tar xzf reference.tar.gz
```

The tarball contains the 36 reference files of the original twelve
cases, the four missing ones among them. It predates the two
`*_mid` cases, which are committed in-tree only. Extracting it on
top of the repository fills in the four missing entries (the other
32 files are byte-identical and harmless to overwrite) and leaves
the `*_mid` files untouched. After
extraction, rerunning `make -C c/tests check` exercises every
comparison.

A companion `reference_pair.tar.gz` is available at
`https://github.com/emol-project/hexpr/releases/download/v1.0.0/reference_pair.tar.gz`
and contains the full `reference_pair/` set. The
`reference_pair/` data is already bundled with the repository, so
the tarball is provided only as a redundant copy for archival
purposes.

## Note: CSF -> determinant expansion

`hexpr_csf_ndets` and `hexpr_csf_to_dets` expand one CSF into the Slater
determinants that span it, at a chosen Ms. This is the seam to any
determinant-basis code.

The full conventions -- genealogical coupling in the Eq. (5) order,
Condon-Shortley CG phase, alpha-then-beta operator order, determinant
ordering, bit-mask layout, normalization -- are stated in
[`c/include/hexpr.h`](include/hexpr.h) and derived in
[`c/src/csf/csf2det.c`](src/csf/csf2det.c). A caller comparing against
another program must match them, or the coefficients will not agree.
To compare against another CSF code, including one that couples the
spins in a different scheme, see
[`docs/COMPARING_WITH_OTHER_CSF_CODES.md`](../docs/COMPARING_WITH_OTHER_CSF_CODES.md).

Two points are worth repeating here because they decide how the API is
used and tested:

- `hexpr_csf_ndets` is **exact, not an upper bound**, unlike
  `hexpr_pair_max_terms` and friends. Buffers sized from it need no
  slack.
- The **overall sign** of a CSF is a convention, not a measured quantity.
  Determinant counts, coefficient magnitudes, unit norm, `<I|J>` and
  `S^2` are all invariant under flipping the sign of a whole CSF, so none
  of them constrains it. Only `validation/csf2det_vs_pyscf.py`, which
  compares the assembled determinant-basis matrix against PySCF, pins the
  sign convention against something outside hexpr. Tests should assert
  *relative* signs, never an individual coefficient's sign.

The `index` argument follows the `hexpr_csf_steps` convention -- the input
position, in `[0, hexpr_drt_ncsfs())`. So do `hexpr_pair_eval` and friends:
there is one index space across the whole API. For a CSFSET sub-DRT the
library maps input positions onto the canonical walk internally, so you
never supply a `hexpr_drt_walk_csfs()` position.

## Memory model

Each handle type (`hexpr_drt_t`, `hexpr_expr_t`) is opaque and
heap-allocated. Ownership rules:

| Call | Returns | Caller frees? | How |
|------|---------|---------------|-----|
| `hexpr_drt_create(...)` | `hexpr_drt_t *` | yes | `hexpr_drt_destroy` |
| `hexpr_drt_from_csfset(norb, csfs, n_csfs)` | `hexpr_drt_t *` | yes | `hexpr_drt_destroy` |
| `hexpr_expr_build(drt)` | `hexpr_expr_t *` | yes | `hexpr_expr_destroy` |
| `hexpr_drt_csfs(drt)` | `const unsigned char *` | **no** | owned by DRT |
| `hexpr_drt_walk_csfs(drt)` | `const unsigned char *` | **no** | owned by DRT |
| `hexpr_last_error()` | `const char *` | **no** | thread-local static |
| `hexpr_version_string()` | `const char *` | **no** | static literal |

`hexpr_expr_read` and `hexpr_drt_node_info` write into caller-allocated
buffers; no allocation crosses the API boundary.

`hexpr_csf_to_dets` likewise writes into caller-allocated buffers. Size
them from `hexpr_csf_ndets`, which is exact rather than an upper bound:
`out_alpha` and `out_beta` each need `ndets * nwords` elements and
`out_coef` needs `ndets`, where `nwords = (norb + 63) / 64`.

Thread safety: read-only queries on a handle are reentrant across threads.
Constructors and destructors are not reentrant on the same handle.
`hexpr_init`/`hexpr_shutdown` must be called from a single thread before
workers start. `hexpr_last_error` is thread-local (C11 `_Thread_local`).

The public API header [`c/include/hexpr.h`](include/hexpr.h) is the
authoritative specification of the API contract.

## Note: internal sub-DRT API

Per-pair / block / subspace evaluation (`hexpr_pair_eval` and friends)
internally builds a small sub-DRT from the `{bra, ket}` CSF set and feeds
it to the same Brooks machinery used for the full DRT. The user-facing
API is unchanged -- `hexpr_drt_create` produces the full DRT and the
pair APIs accept it directly.

For callers that want to construct a sub-DRT directly from a list of
CSF step sequences (e.g., for prototyping, custom analyses, or working
with externally-supplied CSF subsets), the constructor
`hexpr_drt_from_csfset(int norb, const unsigned char *csfs, int n_csfs)`
is available. The returned `hexpr_drt_t *` satisfies all standard DRT
accessors and can be passed to any function that takes a `hexpr_drt_t *`,
including `hexpr_pair_eval`. Validation (step codes in 0..3, common top
node, etc.) is performed on input; failure returns `NULL` and sets
`hexpr_last_error()`.

The kind constant for sub-DRTs is `HEXPR_DRT_CSFSET = 2`
(distinguished from `HEXPR_DRT_SOCI = 0` and `HEXPR_DRT_FOCI = 1`).

**CSF ordering (sub-DRT `csfs`):** for a CSFSET sub-DRT from
`hexpr_drt_from_csfset`, `hexpr_drt_csfs()` / `hexpr_drt_ncsfs()` /
`hexpr_csf_index` / `hexpr_csf_steps` speak the **input set in input
order**: `hexpr_drt_csfs()` returns exactly the CSFs you supplied, in the
order you supplied them, and `hexpr_drt_ncsfs()` is that input count. The
`bra_pos`/`ket_pos` returned by `hexpr_block_eval`/`hexpr_subspace_eval`
index that same input array — consistent with `csfs()`.

The internal DRT may expand to a canonical **superset** (more rows than the
input, when input CSFs share nodes) or dedup (fewer). That canonical view —
in canonical GUGA order (e<u<d<f), the form used internally for evaluation —
is available via **`hexpr_drt_walk_csfs()` / `hexpr_drt_walk_ncsfs()`**.

**CSF ordering and count (SOCI/FOCI):** for a SOCI or FOCI DRT there is no
separate input set, so `hexpr_drt_csfs()` and `hexpr_drt_walk_csfs()` are the
same buffer, and `hexpr_drt_ncsfs()` and `hexpr_drt_walk_ncsfs()` are the same
value. That value is exact: it is the number of complete GUGA walks from the
bottom node to the target state, and `hexpr_drt_create` verifies this count
against its own walk enumeration at construction time -- the two are
guaranteed to agree, or construction fails with `hexpr_last_error()` set
rather than returning a DRT whose count doesn't match its own CSF buffer.
`hexpr_drt_csfs()` (== `hexpr_drt_walk_csfs()`) is allocated to exactly
`hexpr_drt_ncsfs() * hexpr_drt_norb()` bytes and filled that far; sizing a
read of it by `hexpr_drt_ncsfs()` (every binding's `.csfs()` does this) is
safe for every constructible DRT.

For the semantics of `pair_eval` and the per-direction
decomposition convention, see [PAIR_API_SEMANTICS.md](../docs/PAIR_API_SEMANTICS.md).

## See also

- Project root: [`../README.md`](../README.md)
- Public API header: [`c/include/hexpr.h`](include/hexpr.h)
- Index encoding: [`docs/PQRS_INDEX_SPEC.md`](../docs/PQRS_INDEX_SPEC.md)
- Comparing with other CSF codes: [`docs/COMPARING_WITH_OTHER_CSF_CODES.md`](../docs/COMPARING_WITH_OTHER_CSF_CODES.md)
- Python binding: [`../python/README.md`](../python/README.md)
- Fortran binding: [`../fortran/README.md`](../fortran/README.md)
- Julia binding: [`../julia/README.md`](../julia/README.md)
- Ruby binding: [`../ruby/README.md`](../ruby/README.md)
- Methodology paper: [DOI placeholder until preprint is up]
