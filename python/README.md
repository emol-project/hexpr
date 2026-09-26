# hexpr -- Python binding

## Overview

This package provides Python bindings for the **hexpr** C library via the
standard `ctypes` module. No compiled extension module is required -- the
binding loads `libhexpr.so` at runtime, so no separate build step is needed
after installing the library itself.

> macOS: the binding loads `libhexpr.dylib`, and `DYLD_LIBRARY_PATH` replaces
> `LD_LIBRARY_PATH` in the examples below. Resolution is automatic.

The binding exposes the core public API: constructing a Distinct Row Table
(DRT) from either FOCI/SOCI parameters or an explicit CSF set, building
the full energy expression, evaluating individual CSF pairs or rectangular
blocks or square subspaces, and reading `(ij, pqrs, coef)` arrays for
direct use in a CI driver. NumPy arrays are returned for numerical data.

For the broader project context, see the [project root README](../README.md).

## Requirements

- Python 3.9 or later
- `libhexpr.so` -- built from `c/` (see [c/README.md](../c/README.md))
- numpy -- array return types for CSF list and expression terms

## Installation

Build the C library first, then install the Python package in editable mode:

```sh
# from the repository root
make -C c
pip install -e python/
```

`libhexpr.so` must be discoverable at import time. The package searches in
this order:

1. `HEXPR_LIBRARY` environment variable (full path to the `.so`)
2. `../../c/libhexpr.so` relative to the package directory (in-tree build)
3. System linker path / `LD_LIBRARY_PATH`

The simplest way to use the in-tree build without installing to a system
path:

```sh
LD_LIBRARY_PATH=$PWD/c python -c "import hexpr; print(hexpr.version())"
```

## Quick start

```python
import hexpr
import numpy as np

hexpr.init()   # idempotent; safe to call multiple times

# Build a DRT from FOCI/SOCI parameters
drt = hexpr.Drt(
    kind='soci',   # or 'foci'
    nve=4,
    spin_x2=0,     # 2*S: 0=singlet, 2=triplet, ...
    ncore=0,
    nval=4,
    next_=0,       # `next_` avoids collision with Python builtin `next`
)
print(drt.norb, drt.nnodes, drt.ncsfs)

# CSF list as numpy array of step codes (0=E, 1=U, 2=D, 3=F)
csfs = drt.csfs()    # shape (ncsfs, norb), dtype uint8

# Alternatively: build a (sub-)DRT directly from an explicit CSF set.
# See "CSF step codes" below for the meaning of each code.
# All CSFs must reach the same top node -- i.e. the same total nve and
# the same spin coupling at the rightmost orbital.
csf_seq = np.array([
    [3, 3, 0],   # F F E
    [3, 1, 2],   # F U D
    [3, 0, 3],   # F E F
], dtype=np.uint8)
sub = hexpr.Drt.from_csfset(csf_seq, norb=3)
print(sub.norb, sub.ncsfs, sub.nnodes)

# Build the full energy expression (all CSF pairs at once)
expr = hexpr.Expression(drt)
print(expr.nterms_full(), expr.nterms_one())

# Bulk read -- all indices are 0-based
ij, pqrs, coef = expr.terms_full()   # shape (nterms_full,) each
ij1, pqrs1, coef1 = expr.terms_one()

# Or: evaluate selected CSF pairs / blocks / subspaces on demand.
# Low-level (DRT-indexed) workflow:
ij_p, pqrs_p, coef_p = drt.pair_eval(bra=0, ket=1, which='full')
ij_b, pqrs_b, coef_b, bra_pos, ket_pos = drt.block_eval(
    bra_indices=[0, 1, 2], ket_indices=[3, 4, 5], which='full')
ij_s, pqrs_s, coef_s, bra_pos_s, ket_pos_s = drt.subspace_eval(
    indices=[0, 1, 2], which='full')

# High-level (step-code-indexed) workflow:
# Just pass two CSF step-code arrays; the binding builds a transient
# sub-DRT internally and discards it.
bra_steps = np.array([3, 3, 0, 0], dtype=np.uint8)   # F F E E
ket_steps = np.array([3, 0, 3, 0], dtype=np.uint8)   # F E F E
ij, pqrs, coef = hexpr.csf_pair_eval(bra_steps, ket_steps, which='full')
```

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
when constructing a sub-DRT via `Drt.from_csfset` or when calling any of
the high-level `csf_*_eval` functions.

## Index base

`ij` and `pqrs` returned by `terms_full()` / `terms_one()` and by the
pair / block / subspace evaluators are **0-based**, matching the binary
API (`hexpr_expr_read`).

`ij` is a 0-based **packed triangle index** of the `(bra, ket)` CSF pair,
not a row index. Decode it before use:

```python
I = int((np.sqrt(8 * ij + 1) - 1) // 2)   # 0-based bra, I >= J
J = ij - I * (I + 1) // 2                  # 0-based ket
```

then add 1 in-language if your array framework is 1-based. This 0-based
convention is exactly what the PySCF / Psi4 / Fermi.jl example seams
decode (they `+1` on the decoded `(I, J)` for their own indexing).

The text dump format (`reference/*.expr.txt`) adds 1 to `ij` for
compatibility with the Ruby reference output, but does **not** adjust
`pqrs` (a deliberate dump-only `ij`-1 / `pqrs`-0 asymmetry).

For full `(ij, pqrs)` decode logic, see
[`docs/PQRS_INDEX_SPEC.md`](../docs/PQRS_INDEX_SPEC.md).

## API reference

### Module-level

| Name | Kind | Description |
|------|------|-------------|
| `hexpr.init()` | function | Initialize the library (idempotent) |
| `hexpr.shutdown()` | function | Release library resources |
| `hexpr.version()` | function | Return version string |
| `hexpr.Drt(kind, nve, spin_x2, ncore, nval, next_)` | class | Build a DRT from FOCI/SOCI parameters |
| `hexpr.Drt.from_csfset(csfs, norb)` | classmethod | Build a (sub-)DRT from an explicit CSF set |
| `hexpr.Expression(drt)` | class | Build the full energy expression |
| `hexpr.csf_pair_eval(bra_steps, ket_steps, which='full')` | function | High-level pair evaluation from step-code sequences |
| `hexpr.csf_block_eval(bra_csfs, ket_csfs, which='full')` | function | High-level rectangular block evaluation |
| `hexpr.csf_subspace_eval(csfs, which='full')` | function | High-level square subspace evaluation |
| `hexpr.step_parse(s, norb)` | function | Parse `"eudf..."` string to a `uint8` step-code array of length `norb` |

### Drt properties and methods

| Name | Kind | Description |
|------|------|-------------|
| `drt.norb` / `.nnodes` / `.ncsfs` | property | DRT scalar queries |
| `drt.ncore` / `.nval` / `.next_` / `.nve_state` / `.spin_x2` | property | DRT scalar queries (`int`) |
| `drt.kind` | property | DRT kind string (`'soci'`, `'foci'`, or `'csfset'`) |
| `drt.csfs()` | method | CSF step-code array, shape `(ncsfs, norb)`, dtype `uint8` (CSFSET: input order) |
| `drt.walk_csfs()` | method | Canonical GUGA-walk CSF array, shape `(walk_ncsfs, norb)`, `uint8` (CSFSET: the superset/canonical view) |
| `drt.walk_ncsfs` | property | Canonical CSF count (`int`); CSFSET: may differ from `ncsfs` (superset M≥N or dedup) |
| `drt.csf_index(steps)` | method | CSF row index for a given step-code sequence |
| `drt.csf_steps(index)` | method | CSF step codes for a row index; `uint8` array of length `norb` |
| `drt.csf_ndets(index, ms_x2)` | method | Number of determinants CSF `index` expands into at `ms_x2`; exact, not a bound |
| `drt.csf_to_dets(index, ms_x2)` | method | Determinant expansion; returns `(alpha, beta, coef)` -- see [CSF to determinants](#csf-to-determinants) |
| `drt.node_info(index)` | method | Node details `(orb, nve, is)` by index |
| `drt.pair_count(bra, ket, which='full')` | method | Term count for one pair |
| `drt.pair_max_terms(which='full')` | method | Upper bound on terms for any single pair |
| `drt.pair_eval(bra, ket, which='full')` | method | Pair evaluation; returns `(ij, pqrs, coef)` |
| `drt.block_count(bra_indices, ket_indices, which='full')` | method | Term count for a rectangular block |
| `drt.block_max_terms(n_bra, n_ket, which='full')` | method | Upper bound on terms for a block |
| `drt.block_eval(bra_indices, ket_indices, which='full')` | method | Block evaluation; returns `(ij, pqrs, coef, bra_pos, ket_pos)` |
| `drt.subspace_count(indices, which='full')` | method | Term count for a square subspace |
| `drt.subspace_max_terms(n, which='full')` | method | Upper bound on terms for a subspace |
| `drt.subspace_eval(indices, which='full')` | method | Subspace evaluation; returns `(ij, pqrs, coef, bra_pos, ket_pos)` |
| `drt.show(file=None)` | method | Print the DRT to stdout (`file=None`); raises `NotImplementedError` if `file` is given |

### Expression methods

| Name | Kind | Description |
|------|------|-------------|
| `expr.nterms_full()` / `.nterms_one()` | method | Term counts |
| `expr.terms_full()` | method | Full expression: `(ij, pqrs, coef)` arrays |
| `expr.terms_one()` | method | 1-electron subset: `(ij, pqrs, coef)` arrays |
| `expr.show_full(file=None)` | method | Print the full expression to stdout |
| `expr.show_one(file=None)` | method | Print the 1-electron expression to stdout |

For full parameter and return-type details, consult the Python docstrings
(`help(hexpr.Drt)`, etc.) and the C header
[`c/include/hexpr.h`](../c/include/hexpr.h).

### Notes on the block / subspace return values

`block_eval` and `subspace_eval` return five arrays. The extra `bra_pos`
and `ket_pos` arrays give, for each output term, the **position within
the caller's input index arrays** (0-based), not the raw DRT CSF index.
For example, if you call `drt.block_eval(bra_indices=[5, 9, 12],
ket_indices=[3, 4])`, then `bra_pos[k]` is in `{0, 1, 2}` (referring to
the position in `[5, 9, 12]`), and `ket_pos[k]` is in `{0, 1}`. To
recover the raw CSF index, the caller can do
`bra_indices[bra_pos[k]]` and `ket_indices[ket_pos[k]]`.

`csf_block_eval` and `csf_subspace_eval` (the high-level variants) use
the same convention: `bra_pos[k]` and `ket_pos[k]` refer to row indices
into the input `bra_csfs` / `ket_csfs` / `csfs` matrices.

## CSF to determinants

A CSF is a fixed linear combination of Slater determinants. `csf_ndets`
reports how many determinants a CSF expands into at a chosen `Ms`, and
`csf_to_dets` returns those determinants with their coefficients. This is
the seam to any determinant-basis code: it lets a CSF-basis quantity be
re-expressed in the determinant basis, or a determinant-basis wavefunction
be projected onto CSFs.

```python
drt = hexpr.Drt(kind='soci', nve=2, spin_x2=0, ncore=1, nval=1, next_=2)
i = drt.csf_index(hexpr.step_parse("udee", drt.norb))

n = drt.csf_ndets(i, 0)               # ms_x2 = 0; here 2
alpha, beta, coef = drt.csf_to_dets(i, 0)
# alpha.shape == beta.shape == (n, nwords), dtype uint64
# coef.shape  == (n,),                      dtype float64

p = 0                                  # spatial orbital, 0-based
has_alpha = (alpha[0, p >> 6] >> (p & 63)) & 1
```

`ms_x2` is `2*Ms`. Any `|Ms| <= S` with the same parity as `spin_x2` is
accepted, negative values included -- `Ms` is not assumed to equal `S`.

`csf_ndets` is **exact, not an upper bound**: it is
`C(n_open, n_alpha_open)`, so `csf_to_dets` writes exactly that many
determinants and buffers can be sized from it with no slack. This is
unlike `pair_max_terms` and friends, which are bounds.

The `index` argument follows the same convention as `csf_steps`: the input
position, in `0 .. ncsfs-1`. So do `pair_eval` / `block_eval` /
`subspace_eval` -- there is one index space across the whole API. For a
CSFSET sub-DRT the library maps input positions onto the canonical walk
internally, so you never supply a `walk_csfs` position.

### Determinant encoding

`alpha[k]` and `beta[k]` are occupation bit masks for determinant `k`,
split into 64-bit words: bit `p % 64` of word `p // 64` is spatial orbital
`p`, with `nwords = (norb + 63) // 64`. The second axis is kept even when
`nwords == 1`, because collapsing it would change the caller's indexing at
`norb > 64`.

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

All five bindings expose the same `(ndets, nwords)` shape and the same
conventions: the **data shape follows the C contract** rather than each
language's idiom. Only error reporting is adapted per language --
exceptions here, `error()` in Julia, `raise` in Ruby, a status return in
Fortran and C.

## Note: internal sub-DRT mechanism

The `Drt.from_csfset(csfs, norb)` classmethod wraps the C function
`hexpr_drt_from_csfset`, which constructs a DRT whose nodes and arcs are
restricted to those reachable from the supplied CSF set. This is the same
internal mechanism used by the C library's per-pair / block / subspace
evaluation APIs (`hexpr_pair_eval` and friends), exposed here for callers
who want to construct a sub-DRT directly for prototyping, custom
analyses, or working with externally-supplied CSF subsets.

The high-level `hexpr.csf_pair_eval` / `csf_block_eval` /
`csf_subspace_eval` functions wrap this pattern in one step: they accept
step-code sequences, build a transient sub-DRT internally, evaluate, and
discard the sub-DRT. Use them when you do not need to reuse a Drt across
many evaluations.

The returned `Drt` instance from `Drt.from_csfset` satisfies all standard
DRT accessors (`norb`, `nnodes`, `ncsfs`, `csfs()`, etc.) and can be
passed to `hexpr.Expression` in the usual way. Its `kind` is `'csfset'`.

**CSF ordering (sub-DRT `csfs()`):** for a CSFSET sub-DRT from
`Drt.from_csfset` (or the high-level `csf_*_eval`), `csfs()` / `ncsfs` /
`csf_index` / `csf_steps` speak the **input set in input order**: `csfs()`
returns exactly the CSFs you supplied, in the order you supplied them, and
`ncsfs` is that input count. The `bra_pos`/`ket_pos` returned by
`block_eval`/`subspace_eval` index that same input array — consistent with
`csfs()`.

The internal DRT may expand to a canonical **superset** (more rows than the
input, when input CSFs share nodes) or dedup (fewer). That canonical view —
in canonical GUGA order (e<u<d<f), the form used internally for evaluation —
is available via **`walk_csfs()` / `walk_ncsfs`**.

Validation of the input CSFs (step codes in 0..3, common top node) is
performed on the C side; invalid input raises `hexpr.HexprError` with a
descriptive message.

## Building from source

No separate build step is needed beyond what is described in
[Installation](#installation). `pip install -e python/` reads
`python/pyproject.toml` and installs the package in editable mode; the
binding loads `libhexpr.so` at runtime via ctypes, so there is no C
extension to compile.

## Running tests

```sh
# from python/ -- bypasses any PYTHONPATH pollution
cd python
PYTHONPATH=$PWD python3 -m pytest tests/ -v
```

Tests verify DRT structure, the full energy expression on all 12
reference cases, the `Drt.from_csfset` constructor, the per-pair /
block / subspace evaluators (including aggregation against the full
expression), the high-level `csf_*_eval` functions against their
low-level counterparts, and the determinant expansion
(`tests/test_csf2det.py`: determinant counts, bit masks, magnitudes and
relative signs, plus a 70-orbital case that exercises the multi-word
path where `nwords > 1`).

## PYTHONPATH caution

If your shell's `PYTHONPATH` contains `.` or `..` (common in some Intel
oneAPI environments), Python's `PathFinder` will find the repo root
directory as a namespace package named `hexpr` before the
editable-install finder can intercept -- making `import hexpr` resolve
to an empty namespace package with no `init` attribute.

Symptom:
AttributeError: module 'hexpr' has no attribute 'init'
$ python3 -c "import hexpr; print(hexpr.path)"
_NamespacePath(['.../hexpr'])   <-- wrong; should be .../python/hexpr

Fix: run pytest from `python/` with `PYTHONPATH=$PWD` as shown above,
or set `HEXPR_LIBRARY` and invoke with a sanitized path:

```sh
PYTHONPATH=/opt/intel/.../pythonapi:/home/you/src/python \
  python3 -m pytest python/tests/ -v
```

## Dependencies

| Package | Purpose |
|---------|---------|
| numpy   | Array return types for CSF list and expression terms |
| pytest  | Test runner (dev dependency only) |

## See also

- Project root: [`../README.md`](../README.md)
- C library: [`../c/README.md`](../c/README.md)
- Public API header: [`../c/include/hexpr.h`](../c/include/hexpr.h)
- Index encoding: [`../docs/PQRS_INDEX_SPEC.md`](../docs/PQRS_INDEX_SPEC.md)
- Fortran binding: [`../fortran/README.md`](../fortran/README.md)
- Julia binding: [`../julia/README.md`](../julia/README.md)
- Ruby binding: [`../ruby/README.md`](../ruby/README.md)
- Methodology paper: [DOI placeholder until preprint is up]
