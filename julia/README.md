# hexpr -- Julia binding

## Overview

This directory provides a Julia binding for the **hexpr** C library via
Julia's built-in `ccall` mechanism. There is no native code to compile
beyond `libhexpr.so` itself -- the module `julia/src/HExpr.jl` calls the
shared library directly with no wrapper layer.

> macOS: the binding loads `libhexpr.dylib` (via `Libdl.dlext`), and
> `DYLD_LIBRARY_PATH` replaces `LD_LIBRARY_PATH`. Resolution is automatic.

For the broader project context, see the [project root README](../README.md).

## Requirements

- Julia 1.12 or later
- `libhexpr.so` -- built from `c/` (see [c/README.md](../c/README.md))
- Standard library: `Test` (tests)

## Installation

The Julia environment is project-local. Run once from the repository root:

```sh
cd julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

After that, pass `--project=.` (from `julia/`) or `--project=julia` (from
the repository root) on every `julia` invocation. No further step is needed.

## Quick start

```julia
using HExpr

HExpr.init()

# ... see workflows below ...

HExpr.shutdown()
```

### 1. SOCI / FOCI workflow

Build the DRT for a FOCI or SOCI space and read the full expression in
one shot (this example uses the `s_full_4o4e` reference case: 4 valence
electrons, singlet SOCI):

```julia
drt = HExpr.Drt("soci"; nve=4, spin_x2=0, ncore=0, nval=4, next=0)
println("norb=", HExpr.norb(drt), "  ncsfs=", HExpr.ncsfs(drt))  # 4  20

# Build the expression and read all terms
expr = HExpr.Expression(drt)
ij, pqrs, coef = HExpr.read_full(expr)
println("nterms_full=", length(coef))  # 638
```

Save the snippet (including the `using HExpr` / `init` / `shutdown` from
above) as `quickstart.jl` and run:

```sh
cd julia
julia --project=. quickstart.jl
```

### 2. CSF set workflow (per-pair / block / subspace evaluation)

Pass a chosen CSF set directly and evaluate individual pairs, rectangular
blocks, or square subspaces. Each step code is one of `0=E`, `1=U`,
`2=D`, `3=F` (see "CSF step codes" below). All CSFs must reach the same
top node.

```julia
csfs = UInt8[
    3 3 0 0
    3 1 2 0
    3 0 3 0
]

# Low-level: build the DRT once, then evaluate via row indices
sub = HExpr.Drt(csfs, 4)
ij_p, pqrs_p, coef_p = HExpr.pair_eval(sub, 0, 1; which=:full)
ij_b, pqrs_b, coef_b, bra_pos, ket_pos = HExpr.block_eval(
    sub, [0, 1], [1, 2]; which=:full)
ij_s, pqrs_s, coef_s, bra_pos_s, ket_pos_s = HExpr.subspace_eval(
    sub, [0, 1, 2]; which=:full)

# High-level: skip the DRT entirely and pass step-code arrays directly.
# The binding builds a transient sub-DRT internally and discards it.
bra_steps = UInt8[3, 3, 0, 0]
ket_steps = UInt8[3, 0, 3, 0]
ij, pqrs, coef = HExpr.csf_pair_eval(bra_steps, ket_steps; which=:full)
```

## Index base

`ij` and `pqrs` are **0-based**, matching the C binary API. `ij` is a 0-based
*packed triangle index* of the `(bra, ket)` CSF pair (not a row index) — decode
it before use:

```julia
I = Int(floor((sqrt(8 * ij + 1) - 1) / 2))   # 0-based bra, I >= J
J = ij - I * (I + 1) ÷ 2                       # 0-based ket
```

then add 1 for Julia's 1-based arrays (`H[I+1, J+1]`), as `run_fermi.jl` does.
`pqrs` is 0-based too. The text dump (`reference/*.expr.txt`) adds 1 to `ij`
only (to match the Ruby reference corpus); `pqrs` stays 0-based.

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
when constructing a sub-DRT via `Drt(csfs, norb)` or when calling any of
the high-level `csf_*_eval` functions.

## API reference

The module `HExpr` (source: [`julia/src/HExpr.jl`](src/HExpr.jl)) exports
the following public symbols.

**Constants**

| Name | Value | Notes |
|------|-------|-------|
| `HEXPR_OK` | 0 | Success status |
| `HEXPR_ERR_INVALID_ARG` ... `HEXPR_ERR_INTERNAL` | 1...99 | Error codes |
| `HEXPR_DRT_SOCI` / `HEXPR_DRT_FOCI` / `HEXPR_DRT_CSFSET` | 0 / 1 / 2 | DRT kind |
| `HEXPR_EXPR_FULL` / `HEXPR_EXPR_ONE` | 0 / 1 | Expression selector |
| `HEXPR_STEP_E/U/D/F` | 0 to 3 | CSF step codes |

**Lifecycle**

| Name | Description |
|------|-------------|
| `init()` | Initialize library; errors on failure |
| `shutdown()` | Release library resources |
| `version_string()` | Version as `String` |
| `version_major/minor/patch()` | Version integers (`Cint`) |
| `last_error()` | Last error message as `String` |

**Types**

| Type | Description |
|------|-------------|
| `Drt` | Distinct Row Table handle; owns C pointer, finalized automatically |
| `Expression` | Expression handle (not `Expr` -- see [Caveats](#caveats)); owns C pointer |

`Drt` constructors:
- `Drt(kind::Integer, nve, spin_x2, ncore, nval, next)` -- integer-kind variant
- `Drt(kind::AbstractString; nve, spin_x2, ncore, nval, next)` -- string-kind
  convenience (`"soci"`, `"foci"`, or `"csfset"`)
- `Drt(csfs::AbstractMatrix{UInt8}, norb::Integer)` -- build a (sub-)DRT
  from a directly supplied CSF set. `csfs` is a `(n_csfs, norb)` matrix
  of step codes (see "CSF step codes")

**Drt accessors**

| Name | Returns | Notes |
|------|---------|-------|
| `norb(d)`, `nnodes(d)`, `ncsfs(d)` | `Cint` | Dimensions |
| `kind(d)` | `Cint` | `HEXPR_DRT_SOCI`, `HEXPR_DRT_FOCI`, or `HEXPR_DRT_CSFSET` |
| `ncore(d)`, `nval(d)` | `Cint` | Orbital partition |
| `next_orb(d)` | `Cint` | External orbital count (not `next` -- see [Caveats](#caveats)) |
| `nve_state(d)`, `spin_x2(d)` | `Cint` | Electron count, 2*spin |
| `csfs(d)` | `Matrix{UInt8}` | CSF step-code matrix of shape `(ncsfs, norb)` (CSFSET: input order) |
| `walk_csfs(d)` | `Matrix{UInt8}` | Canonical GUGA-walk CSF matrix `(walk_ncsfs, norb)` (CSFSET: the superset/canonical view) |
| `walk_ncsfs(d)` | `Cint` | Canonical CSF count; CSFSET: may differ from `ncsfs` (superset M≥N or dedup) |
| `csf_index(d, steps)` | `Cint` | 0-based CSF row index for a given step-code vector |
| `csf_steps(d, index)` | `Vector{UInt8}` | CSF step codes for a row index; length `norb` |
| `csf_ndets(d, index, ms_x2)` | `Int` | Determinants CSF `index` expands into at `ms_x2`; exact, not a bound |
| `csf_to_dets(d, index, ms_x2)` | `(alpha, beta, coef)` | Determinant expansion -- see [CSF to determinants](#csf-to-determinants) |
| `step_parse(s, norb)` | `Vector{UInt8}` | Parse `"eudf..."` string to step codes (standalone util; placed here next to `csf_steps`) |
| `drt_node_info(d, node_index)` | `NamedTuple` | Node details `(orb, nve, is)` by index |
| `drt_show(d)` | `nothing` | Print the DRT to stdout |

**Expression accessors**

| Name | Returns | Description |
|------|---------|-------------|
| `nterms_full(e)` | `Int64` | Total term count |
| `nterms_one(e)` | `Int64` | 1-electron term count |
| `read_full(e)` | `(ij, pqrs, coef)` | Bulk read -- all terms; `Vector{Int32}`, `Vector{Float64}` |
| `read_one(e)` | `(ij, pqrs, coef)` | Bulk read -- 1-electron subset |
| `expr_show_full(e)` | `nothing` | Print the full expression to stdout |
| `expr_show_one(e)` | `nothing` | Print the 1-electron expression to stdout |

**Per-pair / block / subspace evaluators (low level)**

These take a `Drt` and CSF row indices (0-based). `which` is a Symbol:
`:full` or `:one`, defaulting to `:full`.

| Name | Returns | Description |
|------|---------|-------------|
| `pair_count(d, bra, ket; which=:full)` | `Int` | Term count for one pair |
| `pair_max_terms(d; which=:full)` | `Int` | Upper bound on terms for any single pair |
| `pair_eval(d, bra, ket; which=:full)` | `(ij, pqrs, coef)` | Pair evaluation |
| `block_count(d, bra_idx, ket_idx; which=:full)` | `Int` | Term count for a rectangular block |
| `block_max_terms(d, n_bra, n_ket; which=:full)` | `Int` | Upper bound on terms for a block |
| `block_eval(d, bra_idx, ket_idx; which=:full)` | `(ij, pqrs, coef, bra_pos, ket_pos)` | Block evaluation |
| `subspace_count(d, idx; which=:full)` | `Int` | Term count for a square subspace |
| `subspace_max_terms(d, n; which=:full)` | `Int` | Upper bound on terms for a subspace |
| `subspace_eval(d, idx; which=:full)` | `(ij, pqrs, coef, bra_pos, ket_pos)` | Subspace evaluation |

**High-level evaluators (DRT hidden)**

These take CSF step-code arrays directly and build a transient sub-DRT
internally.

| Name | Returns | Description |
|------|---------|-------------|
| `csf_pair_eval(bra_steps, ket_steps; which=:full)` | `(ij, pqrs, coef)` | Pair evaluation from step-code vectors |
| `csf_block_eval(bra_csfs, ket_csfs; which=:full)` | `(ij, pqrs, coef, bra_pos, ket_pos)` | Block evaluation from step-code matrices |
| `csf_subspace_eval(csfs; which=:full)` | `(ij, pqrs, coef, bra_pos, ket_pos)` | Subspace evaluation from a step-code matrix |

### Notes on the block / subspace return values

`block_eval` and `subspace_eval` return five arrays. The extra `bra_pos`
and `ket_pos` arrays give, for each output term, the **position within
the caller's input index arrays** (0-based), not the raw DRT CSF index.
For example, if you call `block_eval(d, [5, 9, 12], [3, 4])`, then
`bra_pos[k]` is in `{0, 1, 2}` (referring to the position in `[5, 9, 12]`)
and `ket_pos[k]` is in `{0, 1}`. To recover the raw CSF index, the
caller can do `bra_idx[bra_pos[k] + 1]` and `ket_idx[ket_pos[k] + 1]`
(adding 1 for Julia's 1-based indexing).

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

```julia
drt = HExpr.Drt("soci"; nve=2, spin_x2=0, ncore=1, nval=1, next=2)
i = HExpr.csf_index(drt, HExpr.step_parse("udee", HExpr.norb(drt)))

n = HExpr.csf_ndets(drt, i, 0)          # ms_x2 = 0; here 2
alpha, beta, coef = HExpr.csf_to_dets(drt, i, 0)
# size(alpha) == size(beta) == (n, nwords),  Matrix{UInt64}
# length(coef) == n,                         Vector{Float64}

p = 0                                    # spatial orbital, 0-based
has_alpha = (alpha[1, (p >> 6) + 1] >> (p & 63)) & 0x01
```

Note the `+ 1` on the word index and on the determinant index: the
orbital number `p` is 0-based (as every index in this binding is), while
the array subscripts are Julia's 1-based.

`ms_x2` is `2*Ms`. Any `|Ms| <= S` with the same parity as `spin_x2` is
accepted, negative values included -- `Ms` is not assumed to equal `S`.

`csf_ndets` is **exact, not an upper bound**: it is
`C(n_open, n_alpha_open)`, so `csf_to_dets` returns exactly that many
determinants. This is unlike `pair_max_terms` and friends, which are
bounds.

The `index` argument follows the same convention as `csf_steps`: the input
position, in `0 .. ncsfs-1`. So do `pair_eval` / `block_eval` /
`subspace_eval` -- there is one index space across the whole API. For a
CSFSET sub-DRT the library maps input positions onto the canonical walk
internally, so you never supply a `walk_csfs` position.

### Determinant encoding

`alpha[k, :]` and `beta[k, :]` are occupation bit masks for determinant
`k`, split into 64-bit words: bit `p % 64` of word `p ÷ 64` is spatial
orbital `p`, with `nwords = (norb + 63) >> 6`. One row is one
determinant, matching `csfs(d)` returning one row per CSF. The second
axis is kept even when `nwords == 1`, because collapsing it would change
the caller's indexing at `norb > 64`.

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
- **Normalization**: `sum(coef .^ 2) == 1` for any one CSF at every `Ms`.
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
`error()` here, exceptions in Python, `raise` in Ruby, a status return in
Fortran and C.

## Note: internal sub-DRT mechanism

The `Drt(csfs, norb)` constructor wraps the C function
`hexpr_drt_from_csfset`, which builds a DRT whose nodes and arcs are
restricted to those reachable from the supplied CSF set. This is the same
internal mechanism used by the C library's per-pair / block / subspace
evaluation APIs (`hexpr_pair_eval` and friends), exposed here for callers
who want to construct a sub-DRT directly for prototyping, custom analyses,
or working with externally-supplied CSF subsets.

The high-level `csf_pair_eval` / `csf_block_eval` / `csf_subspace_eval`
functions wrap this pattern in one step: they accept step-code arrays,
build a transient sub-DRT internally, evaluate, and discard the sub-DRT.
Use them when you do not need to reuse a Drt across many evaluations.

The returned `Drt` instance satisfies all standard DRT accessors and can
be passed to `HExpr.Expression` in the usual way. Its `kind` is
`HEXPR_DRT_CSFSET`.

**CSF ordering (sub-DRT `csfs`):** for a CSFSET sub-DRT from
`Drt(csfs, norb)` (or the high-level `csf_*_eval`), `csfs(d)` / `ncsfs(d)`
/ `csf_index` / `csf_steps` speak the **input set in input order**:
`csfs(d)` returns exactly the CSFs you supplied, in the order you supplied
them, and `ncsfs(d)` is that input count. The `bra_pos`/`ket_pos` returned
by `block_eval`/`subspace_eval` index that same input array — consistent
with `csfs(d)`.

The internal DRT may expand to a canonical **superset** (more rows than the
input, when input CSFs share nodes) or dedup (fewer). That canonical view —
in canonical GUGA order (e<u<d<f), the form used internally for evaluation —
is available via **`walk_csfs(d)` / `walk_ncsfs(d)`**.

Validation of the input CSFs (step codes in 0..3, common top node) is
performed on the C side; invalid input raises an error with a descriptive
message from `last_error()`.

## Building from source

There is no separate build step. `Pkg.instantiate()` (see
[Installation](#installation)) downloads and precompiles the registered
dependencies. The `HExpr` module is pure Julia.

## Running tests

```sh
julia --project=julia -e 'using Pkg; Pkg.test()'
# Or from julia/ directory:
# cd julia && julia --project=. -e 'using Pkg; Pkg.test()'
```

Tests pass across six test sets (`smoke`, `drt`, `expr`, `csfset`,
`pair`, `completeness`) covering lifecycle, DRT construction and
accessors, expression build and read on the `s_full_4o4e` reference case,
the `Drt(csfs, norb)` constructor, the per-pair / block / subspace
evaluators (including aggregation against the full expression), the
high-level `csf_*_eval` functions against their low-level counterparts,
and the binding-completeness wrappers. The determinant expansion is
covered in `test/test_completeness.jl`: determinant counts, bit masks,
magnitudes and relative signs, plus a 70-orbital case that exercises the
multi-word path where `nwords > 1`.

## See also

- Project root: [`../README.md`](../README.md)
- C library: [`../c/README.md`](../c/README.md)
- Public API header: [`../c/include/hexpr.h`](../c/include/hexpr.h)
- Index encoding: [`../docs/PQRS_INDEX_SPEC.md`](../docs/PQRS_INDEX_SPEC.md)
- Python binding: [`../python/README.md`](../python/README.md)
- Fortran binding: [`../fortran/README.md`](../fortran/README.md)
- Ruby binding: [`../ruby/README.md`](../ruby/README.md)
- Methodology paper: [DOI placeholder until preprint is up]
- Full examples: [`../example_ci/`](../example_ci/) -- end-to-end CI
  calculations in Julia (H2 / STO-3G CAS(2,2) and H2O / 6-31G(d) SOCI)

## Caveats

**`Expression`, not `Expr`**
Julia's `Base` module exports `Expr` (the AST node type). The expression
handle type is therefore named `Expression` to avoid the clash. Exporting a
name `Expr` from `HExpr` would silently shadow `Base.Expr` for any module
that does `using HExpr` -- naming it `Expression` prevents the conflict.

**`next_orb`, not `next`**
`Base.next` was the iteration protocol in Julia 0.x (deprecated in 1.x,
removed in later releases). The DRT accessor for the external orbital count
is named `next_orb` to avoid any residual conflict with the `Base`
namespace.

**HDF5.jl axis ordering (C row-major vs Julia column-major)**
The example `example_ci/h2o_soci/run.jl` reads arrays written by h5py
(Python / row-major). Julia's HDF5.jl automatically transposes array axes
on read to match Julia's column-major convention. For the symmetric
tensors in that example (`h1eff` and `eri_act` with 8-fold permutation
symmetry) the transposition produces identical values, so no `permutedims`
is needed. Users loading non-symmetric arrays from h5py-written HDF5
files should account for this axis reversal explicitly.

**Floating-point coefficients: use `isapprox`, not `==`**
The first coefficient in the `s_full_4o4e` full expression is
`nextfloat(2.0)` = `2.0000000000000004`, not exactly `2.0`. The deviation
is a single ULP from the integer accumulation path inside the C library.
The test suite uses `isapprox(coef[1], 2.0)` (the `isapprox` function) rather
than `==` to allow for this.
