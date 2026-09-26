# hexpr -- Ruby binding

## Overview

This directory provides a Ruby binding for the **hexpr** C library via the
[`ffi`](https://github.com/ffi/ffi) gem. No native Ruby extension needs to be
compiled -- the `ffi` gem handles C interop at runtime, loading `libhexpr.so`
directly.

> macOS: the binding loads `libhexpr.dylib` (via `FFI::Platform::LIBSUFFIX`),
> and `DYLD_LIBRARY_PATH` replaces `LD_LIBRARY_PATH`. Resolution is automatic.

## Requirements

- Ruby 3.1 or later (tested on 3.1.2)
- `libhexpr.so` -- built from `c/` (see [c/README.md](../c/README.md))
- C compiler and Ruby development headers -- needed only for the `ffi` gem's
  native extension build during `bundle install`, not for hexpr usage itself
- Bundler 2.3 or later

## Installation

```sh
cd ruby
bundle config set --local path 'vendor/bundle'
bundle install
```

The `bundle config` line isolates gem installs under `ruby/vendor/bundle/`
instead of the user's global gem environment. The `vendor/` directory is
gitignored. This step is needed only once; subsequent invocations just require
`bundle exec`.

## Quick start

```ruby
require "hexpr"

HExpr.init

# ... see workflows below ...

HExpr.shutdown
```

### 1. SOCI / FOCI workflow

Build the DRT for a FOCI or SOCI space and read the full expression in
one shot (this example uses the `s_full_4o4e` reference case: 4 valence
electrons, singlet SOCI):

```ruby
drt = HExpr::Drt.new(kind: "soci", nve: 4, spin_x2: 0, ncore: 0, nval: 4, next: 0)
puts "norb=#{drt.norb}  ncsfs=#{drt.ncsfs}"  # norb=4  ncsfs=20

# Build the expression and read all terms
expr = HExpr::Expression.new(drt)
ij, pqrs, coef = expr.read_full
puts "nterms_full=#{ij.length}"  # nterms_full=638
```

Save the snippet (including the `require` / `init` / `shutdown` from
above) as `quickstart.rb` and run:

```sh
cd ruby
bundle exec ruby -Ilib quickstart.rb
```

### 2. CSF set workflow (per-pair / block / subspace evaluation)

Pass a chosen CSF set directly and evaluate individual pairs, rectangular
blocks, or square subspaces. Each step code is one of `0=E`, `1=U`,
`2=D`, `3=F` (see "CSF step codes" below). All CSFs must reach the same
top node.

```ruby
csfs = [
  [3, 3, 0, 0],
  [3, 1, 2, 0],
  [3, 0, 3, 0],
]

# Low-level: build the DRT once, then evaluate via row indices
sub = HExpr::Drt.from_csfset(csfs, 4)
ij_p, pqrs_p, coef_p = sub.pair_eval(0, 1, which: :full)
ij_b, pqrs_b, coef_b, bra_pos, ket_pos = sub.block_eval(
  [0, 1], [1, 2], which: :full)
ij_s, pqrs_s, coef_s, bra_pos_s, ket_pos_s = sub.subspace_eval(
  [0, 1, 2], which: :full)

# High-level: skip the DRT entirely and pass step-code arrays directly.
# The binding builds a transient sub-DRT internally and discards it.
bra_steps = [3, 3, 0, 0]
ket_steps = [3, 0, 3, 0]
ij, pqrs, coef = HExpr.csf_pair_eval(bra_steps, ket_steps, which: :full)
```

## Index base

`ij` and `pqrs` are **0-based**, matching the C binary API. `ij` is a 0-based
*packed triangle index* of the `(bra, ket)` CSF pair (not a row index) — decode
it before use:

```ruby
i = ((Math.sqrt(8 * ij + 1) - 1) / 2).to_i   # 0-based bra, i >= j
j = ij - i * (i + 1) / 2                       # 0-based ket
```

Ruby arrays are 0-based, so `(i, j)` index directly. `pqrs` is 0-based too. The
text dump (`reference/*.expr.txt`) adds 1 to `ij` only (to match the Ruby
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
when constructing a sub-DRT via `Drt.from_csfset` or when calling any of
the high-level `csf_*_eval` methods.

## API reference

Source: [`ruby/lib/hexpr.rb`](lib/hexpr.rb).

**Constants**

| Name | Value | Notes |
|------|-------|-------|
| `HEXPR_OK` | 0 | Success status |
| `HEXPR_ERR_INVALID_ARG` ... `HEXPR_ERR_INTERNAL` | 1...99 | Error codes |
| `HEXPR_DRT_SOCI` / `HEXPR_DRT_FOCI` / `HEXPR_DRT_CSFSET` | 0 / 1 / 2 | DRT kind |
| `HEXPR_EXPR_FULL` / `HEXPR_EXPR_ONE` | 0 / 1 | Expression selector |
| `HEXPR_STEP_E/U/D/F` | 0 to 3 | CSF step codes |

**Module functions** (`HExpr.*`)

| Name | Description |
|------|-------------|
| `HExpr.init` | Initialize library; raises on failure |
| `HExpr.shutdown` | Release library resources |
| `HExpr.version_string` | Version as `String` |
| `HExpr.version_major/minor/patch` | Version integers |
| `HExpr.last_error` | Last error message as `String` |
| `HExpr.csf_pair_eval(bra_steps, ket_steps, which: :full)` | High-level pair evaluation from step-code arrays |
| `HExpr.csf_block_eval(bra_csfs, ket_csfs, which: :full)` | High-level rectangular block evaluation |
| `HExpr.csf_subspace_eval(csfs, which: :full)` | High-level square subspace evaluation |
| `HExpr.step_parse(str, norb)` | Parse `"eudf..."` string to an `Array<Integer>` of step codes, length `norb` |

**Classes**

| Class | Description |
|-------|-------------|
| `HExpr::Drt` | DRT handle; memory freed via `FFI::AutoPointer` on GC |
| `HExpr::Expression` | Expression handle; memory freed via `FFI::AutoPointer` on GC |

`HExpr::Drt` constructors:

- `HExpr::Drt.new(kind:, nve:, spin_x2:, ncore:, nval:, next:)` -- build a
  DRT from FOCI/SOCI parameters. `kind` is a `String` (`"soci"` or `"foci"`)
  or an `Integer` constant. See [Caveats](#caveats) for the `next:` keyword.
- `HExpr::Drt.from_csfset(csfs, norb)` -- build a (sub-)DRT from a directly
  supplied CSF set. `csfs` is an `Array<Array<Integer>>` (rows are CSFs as
  step codes; see "CSF step codes"); `norb` is the orbital count. The
  inverse direction -- reading a DRT's step codes back out as the same
  `Array<Array<Integer>>` -- is `Drt#csfs`.

**`HExpr::Drt` accessors**

| Name | Returns | Notes |
|------|---------|-------|
| `norb`, `nnodes`, `ncsfs` | `Integer` | Dimensions |
| `kind` | `Integer` | `HEXPR_DRT_SOCI`, `HEXPR_DRT_FOCI`, or `HEXPR_DRT_CSFSET` |
| `ncore`, `nval` | `Integer` | Orbital partition |
| `next_orb` | `Integer` | External orbital count (not `next` -- see [Caveats](#caveats)) |
| `nve_state`, `spin_x2` | `Integer` | Electron count, 2*spin |
| `csfs` | `Array<Array<Integer>>` | CSF step-code matrix, `ncsfs` x `norb`, each element in 0..3. Mirrors the input format of `Drt.from_csfset`, and the output round-trips back into it: for a sub-DRT the rows **are** the input set in input order (the canonical/superset view is `walk_csfs`) |
| `walk_csfs` | `Array<Array<Integer>>` | Canonical GUGA-walk CSF matrix, `walk_ncsfs` x `norb` (CSFSET: the superset/canonical view; SOCI/FOCI: equals `csfs`) |
| `walk_ncsfs` | `Integer` | Canonical CSF count; CSFSET: may differ from `ncsfs` (superset M≥N or dedup) |
| `csf_index(steps)` | `Integer` | 0-based CSF row index for a given step-code array |
| `csf_steps(index)` | `Array<Integer>` | CSF step codes for a row index; length `norb` |
| `csf_ndets(index, ms_x2)` | `Integer` | Determinants CSF `index` expands into at `ms_x2`; exact, not a bound |
| `csf_to_dets(index, ms_x2)` | `[alpha, beta, coef]` | Determinant expansion -- see [CSF to determinants](#csf-to-determinants) |
| `node_info(index)` | `Array<Integer>` | Node details `[orb, nve, is]` by index |
| `show` | `nil` | Print the DRT to stdout |

**`HExpr::Drt` per-pair / block / subspace methods (low level)**

These take CSF row indices (0-based). `which:` is a Symbol: `:full` or
`:one`, defaulting to `:full`.

| Name | Returns | Description |
|------|---------|-------------|
| `pair_count(bra, ket, which: :full)` | `Integer` | Term count for one pair |
| `pair_max_terms(which: :full)` | `Integer` | Upper bound on terms for any single pair |
| `pair_eval(bra, ket, which: :full)` | `[ij, pqrs, coef]` | Pair evaluation |
| `block_count(bra_indices, ket_indices, which: :full)` | `Integer` | Term count for a rectangular block |
| `block_max_terms(n_bra, n_ket, which: :full)` | `Integer` | Upper bound on terms for a block |
| `block_eval(bra_indices, ket_indices, which: :full)` | `[ij, pqrs, coef, bra_pos, ket_pos]` | Block evaluation |
| `subspace_count(indices, which: :full)` | `Integer` | Term count for a square subspace |
| `subspace_max_terms(n, which: :full)` | `Integer` | Upper bound on terms for a subspace |
| `subspace_eval(indices, which: :full)` | `[ij, pqrs, coef, bra_pos, ket_pos]` | Subspace evaluation |

**`HExpr::Expression` methods**

| Name | Returns | Description |
|------|---------|-------------|
| `nterms_full` | `Integer` | Total term count |
| `nterms_one` | `Integer` | 1-electron term count |
| `read_full` | `[ij, pqrs, coef]` | Bulk read -- all terms; `Array<Integer>`, `Array<Float>` |
| `read_one` | `[ij, pqrs, coef]` | Bulk read -- 1-electron subset |
| `show_full` | `nil` | Print the full expression to stdout |
| `show_one` | `nil` | Print the 1-electron expression to stdout |

### Notes on the block / subspace return values

`block_eval` and `subspace_eval` return five arrays. The extra `bra_pos`
and `ket_pos` arrays give, for each output term, the **position within
the caller's input index arrays** (0-based), not the raw DRT CSF index.
For example, if you call `drt.block_eval([5, 9, 12], [3, 4])`, then
`bra_pos[k]` is in `[0, 1, 2]` (referring to the position in
`[5, 9, 12]`) and `ket_pos[k]` is in `[0, 1]`. To recover the raw CSF
index, the caller can do `bra_indices[bra_pos[k]]` and
`ket_indices[ket_pos[k]]`.

`HExpr.csf_block_eval` and `HExpr.csf_subspace_eval` (the high-level
variants) use the same convention: `bra_pos[k]` and `ket_pos[k]` refer
to row indices into the input `bra_csfs` / `ket_csfs` / `csfs` arrays.

## CSF to determinants

A CSF is a fixed linear combination of Slater determinants. `csf_ndets`
reports how many determinants a CSF expands into at a chosen `Ms`, and
`csf_to_dets` returns those determinants with their coefficients. This is
the seam to any determinant-basis code: it lets a CSF-basis quantity be
re-expressed in the determinant basis, or a determinant-basis wavefunction
be projected onto CSFs.

```ruby
drt = HExpr::Drt.new(kind: "soci", nve: 2, spin_x2: 0,
                     ncore: 1, nval: 1, next: 2)
i = drt.csf_index(HExpr.step_parse("udee", drt.norb))

n = drt.csf_ndets(i, 0)                 # ms_x2 = 0; here 2
alpha, beta, coef = drt.csf_to_dets(i, 0)
# alpha and beta: Array<Array<Integer>>, n rows x nwords columns
# coef:           Array<Float>, length n

p = 0                                    # spatial orbital, 0-based
has_alpha = (alpha[0][p / 64] >> (p % 64)) & 1
```

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

### Determinant encoding, and why it is not one Integer

`alpha[k]` and `beta[k]` are occupation bit masks for determinant `k`,
split into 64-bit words: bit `p % 64` of word `p / 64` is spatial orbital
`p`, with `nwords = (norb + 63) / 64`.

Ruby's `Integer` is arbitrary-precision, so a whole mask could have been
returned as a single value and read with `alpha[k][p]`. **That was
considered and deliberately not done.** This is a wrapper: its job is to
present the C contract as the C API states it, and the other four
bindings expose the same `(ndets, nwords)` shape. A Ruby-only packed form
would have to be read back against a different contract from every other
language, for a convenience that costs one division and one modulo. A
packed accessor can be added later without breaking this one; collapsing
this one later could not be undone.

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
`raise` here, exceptions in Python, `error()` in Julia, a status return in
Fortran and C.

## Note: internal sub-DRT mechanism

The `Drt.from_csfset(csfs, norb)` class method wraps the C function
`hexpr_drt_from_csfset`, which builds a DRT whose nodes and arcs are
restricted to those reachable from the supplied CSF set. This is the same
internal mechanism used by the C library's per-pair / block / subspace
evaluation APIs (`hexpr_pair_eval` and friends), exposed here for callers
who want to construct a sub-DRT directly for prototyping, custom analyses,
or working with externally-supplied CSF subsets.

The high-level `HExpr.csf_pair_eval` / `csf_block_eval` /
`csf_subspace_eval` functions wrap this pattern in one step: they accept
step-code arrays, build a transient sub-DRT internally, evaluate, and
discard the sub-DRT. Use them when you do not need to reuse a Drt
across many evaluations.

The returned `HExpr::Drt` instance satisfies all standard DRT accessors
and can be passed to `HExpr::Expression.new` in the usual way. Its `kind`
is `HEXPR_DRT_CSFSET`.

**CSF ordering (sub-DRT `csfs`):** for a CSFSET sub-DRT from
`Drt.from_csfset` (or the high-level `csf_*_eval`), `csfs` / `ncsfs` /
`csf_index` / `csf_steps` speak the **input set in input order**: `csfs`
returns exactly the CSFs you supplied, in the order you supplied them, and
`ncsfs` is that input count. The `bra_pos`/`ket_pos` returned by
`block_eval`/`subspace_eval` index that same input array — consistent with
`csfs`.

The internal DRT may expand to a canonical **superset** (more rows than the
input, when input CSFs share nodes) or dedup (fewer). That canonical view —
in canonical GUGA order (e<u<d<f), the form used internally for evaluation —
is available via **`walk_csfs` / `walk_ncsfs`**.

Validation of the input CSFs (step codes in 0..3, common top node) is
performed on the C side; invalid input raises a `RuntimeError` with a
descriptive message from `HExpr.last_error`. Shape validation (the array
is a 2-D `Array<Array<Integer>>` with each row of length `norb`) is
performed in Ruby and raises `ArgumentError`.

## Building from source

There is no separate build step beyond `bundle install` (see
[Installation](#installation)). `bundle install` compiles the `ffi` gem's
C extension automatically. The `HExpr` module itself is pure Ruby.

## Running tests

```sh
cd ruby && bundle exec ruby -Ilib -Itest \
  -e 'Dir["test/test_*.rb"].each { |f| require_relative f }'
```

Or run test files individually:

```sh
cd ruby
bundle exec ruby -Ilib -Itest test/test_smoke.rb
bundle exec ruby -Ilib -Itest test/test_drt.rb
bundle exec ruby -Ilib -Itest test/test_expr.rb
bundle exec ruby -Ilib -Itest test/test_csfs.rb
bundle exec ruby -Ilib -Itest test/test_csfset.rb
bundle exec ruby -Ilib -Itest test/test_pair.rb
bundle exec ruby -Ilib -Itest test/test_completeness.rb
```

`test/helper.rb` puts `lib/` on `$LOAD_PATH` itself, so `-Ilib -Itest` is
belt-and-braces; `ruby test/test_smoke.rb` from `ruby/` works too, as long
as the `ffi` gem is visible to the Ruby being invoked.

Tests cover lifecycle, DRT construction and accessors, expression build
and read on the `s_full_4o4e` reference case, the `Drt.from_csfset`
class method, the per-pair / block / subspace methods (including
aggregation against the full expression), the high-level
`HExpr.csf_*_eval` functions against their low-level counterparts, and
the binding-completeness wrappers. The determinant expansion is covered
in `test/test_completeness.rb`: determinant counts, bit masks,
magnitudes and relative signs, plus a 70-orbital case that exercises the
multi-word path where `nwords > 1`.

## See also

- Project root: [`../README.md`](../README.md)
- C library: [`../c/README.md`](../c/README.md)
- Public API header: [`../c/include/hexpr.h`](../c/include/hexpr.h)
- Index encoding: [`../docs/PQRS_INDEX_SPEC.md`](../docs/PQRS_INDEX_SPEC.md)
- Python binding: [`../python/README.md`](../python/README.md)
- Fortran binding: [`../fortran/README.md`](../fortran/README.md)
- Julia binding: [`../julia/README.md`](../julia/README.md)
- Methodology paper: [DOI placeholder until preprint is up]

## Caveats

**`next:` is valid as a keyword argument (accessor is `next_orb`)**
`next` is a Ruby reserved word, but Ruby 3.1 accepts it as a keyword argument
name in a method signature. The constructor `Drt.new(... next: 0)` works
correctly; inside the constructor body the value is retrieved via
`binding.local_variable_get(:next)` because `next` cannot appear as a bare
identifier. The accessor method for the external orbital count is named
`next_orb` rather than `next`, since a method named `next` would conflict with
`Object#next` and cannot be defined on an arbitrary class.

**Loaded via `-Ilib`, not packaged as a gem**
`hexpr.rb` is not (yet) published as a RubyGem. Load it by passing `-Ilib`
(or `-I path/to/ruby/lib`) on the Ruby command line, or add the `lib/`
directory to `$LOAD_PATH` in your script. Gem packaging is deferred.

**`FFI::AutoPointer` with class-method release; `@parent` keeps `Drt` alive**
Both `Drt` and `Expression` wrap their C pointers in `FFI::AutoPointer` with
a `self.release` class method so that the underlying `hexpr_drt_destroy` /
`hexpr_expr_destroy` is called automatically on GC. `Expression` also holds
`@parent = drt` to prevent the `Drt` object from being collected while the
`Expression` is still alive.

**`bundle config set --local path 'vendor/bundle'` keeps gems local**
Running `bundle install` without this config line installs gems into the
user's global gem environment. The `bundle config` line scopes the install
to `ruby/vendor/bundle/`, which is gitignored. If you skip this step and
install globally, later `bundle exec` invocations will still work, but the
gems will not be isolated to the project.
