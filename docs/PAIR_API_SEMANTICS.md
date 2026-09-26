# Per-pair API semantics: directional decomposition of matrix elements

> The "Contract" section below is normative for hexpr v1.0.0 and later.
> The rest of this document is explanatory: it records how the behavior
> was measured and why it is the way it is. Where the two appear to
> disagree, the contract governs.

This document describes the semantics of the per-pair / block /
subspace evaluation routines exposed by hexpr. It covers
`hexpr_pair_eval`, `hexpr_pair_count`, and their block and
subspace counterparts (C names), and the equivalent high-level
calls in the Python, Julia, and Ruby bindings.

The directional behavior described here is a well-understood
property of the per-pair evaluator, not a defect. Its origin was
characterized by a controlled measurement (see
`docs/DECISIONS/directional-pair-decomposition.md`); the relevant
findings are stated inline below where they matter. The behavior
is easy to mis-read on first contact, so this note makes the
convention — and what is and is not established about it —
explicit.

## Contract

What a caller may rely on, and what it may not. Everything outside this
section is description, not guarantee.

### Guaranteed

- **G1. One index space.** `I` and `J` are INPUT-ORDER CSF indices in
  `[0, hexpr_drt_ncsfs())`, the same space as `hexpr_drt_csfs`,
  `hexpr_csf_index` and `hexpr_csf_steps`. This holds for every kind of
  DRT, including a CSFSET sub-DRT whose canonical walk reorders, extends
  or dedups the input.
- **G2. `ij` is not an index into that space.** The `ij` in a returned
  term is a canonical pair address (`docs/PQRS_INDEX_SPEC.md`). Decoding
  it does not yield arguments for these calls. Attribute a term to a pair
  with the `bra_pos` / `ket_pos` outputs of the block and subspace forms.
- **G3. Every term returned is a true term.** Each `(ij, pqrs, coef)`
  from `hexpr_pair_eval(drt, I, J, which)` is a term of the decomposition
  of `H[I, J]` restricted by `which`. The call returns a subset of that
  decomposition, never a term that does not belong to it.
- **G4. The union of the two argument orders is complete.** For any pair
  `(I, J)`, the union of the terms from `hexpr_pair_eval(drt, I, J, w)`
  and `hexpr_pair_eval(drt, J, I, w)` is the complete decomposition of
  `H[I, J]` under `w`.
- **G5. Coefficients agree on the overlap.** If a `(ij, pqrs)` key is
  returned by both argument orders, it carries the same coefficient in
  both, so the union in G4 is well defined and involves no double
  counting.
- **G6. Block and subspace are sums of single directions.**
  `hexpr_block_eval(drt, B, K, w)` returns, for every ordered pair
  `(B[i], K[j])`, exactly the terms `hexpr_pair_eval(drt, B[i], K[j], w)`
  would return, tagged with `bra_pos = i` and `ket_pos = j`.
  `hexpr_subspace_eval(drt, S, w)` is `hexpr_block_eval(drt, S, S, w)`.
  See "What block and subspace do, exactly" below for what this implies.

### Not guaranteed

These are observed properties of the current implementation, documented
because they are useful, but they are **not** part of the contract and
may change without a major version bump. Code that must keep working
should not depend on them.

- **N1. The strict-subset shape.** That one argument order's term set is
  always a subset of the other's (rather than merely overlapping) is a
  measurement, taken on the cases in this document. Use the dedup form
  below, not the keep-larger form, if you need the complete
  decomposition and want to be insulated from this.
- **N2. Which order is the larger one.** It depends on the pair, and
  nothing in the API tells you in advance.
- **N3. Term counts and term ordering.** Neither the number of terms in
  a single direction nor the order of terms within an output is
  specified. `hexpr_pair_max_terms` bounds the count; it does not
  predict it.
- **N4. The single-direction behavior of the bindings'
  `csf_pair_eval`.** See "Step-based helpers in the bindings" below.

### Performance

Nothing in this contract is a statement about cost. In particular, the
per-pair routines build and discard a two-CSF sub-DRT internally on every
call, and so do the bindings' step-based helpers. Evaluating many pairs
of the same space one at a time therefore repeats that construction; the
block and subspace forms amortize it across the block. That is current
behavior, not a guarantee in either direction.

## Which index space

`I` and `J` below are INPUT-ORDER CSF indices, in
`[0, hexpr_drt_ncsfs())` -- the same space as `hexpr_csf_index` and
`hexpr_csf_steps`. There is one index space across the whole API. For
a CSFSET sub-DRT the library maps them onto the canonical walk
internally; you never supply a `hexpr_drt_walk_csfs()` position.

The `ij` value in the returned terms is the exception: it stays a
canonical pair address (`docs/PQRS_INDEX_SPEC.md`). Use the
`bra_pos` / `ket_pos` outputs of the block and subspace forms to
attribute a term to a pair, and do not decode `ij` and feed the result
back to these calls.

## What hexpr returns for a single pair

For a single off-diagonal CSF pair `(I, J)` with `I != J`, the
Hamiltonian matrix element

    H[I, J] = sum_k coef_k * integral(pqrs_k)

is symmetric in `I` and `J`: `H[I, J] = H[J, I]`. The sum on the
right, however, is not the only way to write down a decomposition.
hexpr emits *one specific decomposition*, selected by the
direction implied by the (bra, ket) ordering of the call.

That direction matters. Concretely:

- `hexpr_pair_eval(drt, I, J, "full")` returns the terms selected
  for the direction in which `I` is the bra.
- `hexpr_pair_eval(drt, J, I, "full")` returns the terms selected
  for the opposite direction.

Which direction you get is decided by the argument order alone.
(Internally the filter compares path weights inside the two-CSF
sub-DRT built for the pair, not the indices you passed; see below.)

The two sets are different decompositions of the same matrix
element. Their term counts can differ, and which direction is
"complete" depends on the pair.

## Empirical structure of the two decompositions

Everything in this section is measurement. G4 and G5 above are the parts
that were promoted to the contract; property 1 below is N1 and was not.

On the s_full_4o4e test case (20 singlet CSFs from a CAS(4,4)
Full-CI DRT), iterating all 190 distinct off-diagonal pairs:

| relationship between the two directions | count |
|---|---|
| forward direction's set is a strict subset of backward's | 137 |
| backward direction's set is a strict subset of forward's | 26 |
| forward and backward have no terms (both zero)           | 27 |
| forward and backward are identical (size > 0)            | 0 |
| neither direction's set is a subset of the other         | 0 |
| overlapping `(ij, pqrs)` key with disagreeing coefficient | 0 |

Two observations follow:

1. **Either direction is a subset of the other.** There is never
   a term that lives only in one direction while the other has
   its own unique terms; one of the two directions always
   contains the full set, and the other contains a subset of it.
   Which direction is full depends on the pair.

2. **Coefficients agree on the overlap.** When a `(ij, pqrs)`
   key appears in both directions, it carries the same
   coefficient. There is no double counting to worry about.

These properties make merging the two decompositions trivial.

## What produces the directional split (measured)

The split is produced by the **direction filter** in
`add_expr_pair` (`c/src/csf/csf_pair_core.c`): `delta = bra - ket`,
then loop elements with `weight_g - weight_f != delta` are
rejected — the comment there calls this rejecting "the
directionally-opposite half of every loop element." This filter,
together with the `bra`-dependent path-weight selection
(`target_sum = bra - weight_g`) and coefficient computation it
gates, is what makes the two orderings yield different sets.

This was confirmed by a controlled counterfactual on s_full_4o4e
(190 pairs):

- **Disabling the direction filter** does not symmetrize the
  output; it destroys the clean table outright — 189 of 190 pairs
  become mutually *incomparable* and coefficients on shared keys
  stop agreeing (differences up to 3.0). The filter is therefore
  load-bearing; direction is carried by the filter plus the
  `bra`-dependent path-weight and coefficient formulas downstream
  of it, not by the filter alone.

- **The strict-subset shape** (property 1 above; `incomparable =
  0`) is a property of this direction filter. It is fully present
  with the internal loop-closure deduplicator (`rhash` in
  `c/src/loopgen/expr_tree.c`) disabled: with `rhash` off, the
  `forward ⊂ backward` count is unchanged (137) and `incomparable`
  stays 0. `rhash` is a minor loop-closure deduplicator whose only
  measured effect on this case is folding a single orientation
  pair (one pair that is otherwise `equal` becomes a strict subset
  under the fold). It is **not** the cause of the subset structure.

- **The no-double-counting property** (property 2) likewise
  survives `rhash` removal and is broken only by removing the
  direction filter, so it too follows from the filter.

What remains an **analytic** statement, not something the above
build-level experiment can prove: *why* one walk direction always
contains the other (a strict subset) rather than merely
overlapping. The measurements establish that the property is
governed by the direction filter and is independent of `rhash`;
the reason it takes the strict-subset form is a property of the
GUGA loop algebra underlying the reduction (see the methodology
paper's reduction / telescope sections), not a coincidence of the
implementation.

## What block and subspace do, exactly

G6 says the block forms are sums of single directions, one per ordered
pair. That has a consequence worth stating plainly, because it is easy
to assume the opposite:

- **`hexpr_subspace_eval(drt, S, w)` is complete for every pair inside
  `S`.** Because bra and ket range over the same set, both `(I, J)` and
  `(J, I)` are evaluated, so G4's union is present in the output. This
  is the form to use for a CI Hamiltonian over a CSF space.
- **`hexpr_block_eval(drt, B, K, w)` is complete for a pair only if that
  pair appears in both orders**, i.e. only for pairs with both members in
  `B` and both in `K`. For a rectangular block with `B` and `K` disjoint
  -- the usual case, an off-diagonal block of a partitioned matrix --
  **every pair is evaluated in one direction only**, and the result is a
  partial decomposition of each element.

For a disjoint off-diagonal block, either call
`hexpr_block_eval(drt, B, K, w)` and `hexpr_block_eval(drt, K, B, w)` and
merge by G5, or evaluate the pairs individually with the helper below.

This distinction does not affect the summed matrix element when the
caller symmetrizes correctly over a full space, which is why it went
unnoticed: `hexpr_subspace_eval` over the whole space is complete, and
that is what the reference comparisons exercise.

## Step-based helpers in the bindings

The Python, Julia and Ruby bindings expose `csf_pair_eval(bra_steps,
ket_steps, which)`, which takes two step-code sequences instead of
indices. It builds a two-CSF sub-DRT from the two rows, calls
`hexpr_pair_eval` once, and discards the sub-DRT. There is no C-level
equivalent.

Two things about it are worth knowing:

- **It returns one direction, not the complete decomposition.** The name
  does not say so. It is the `(bra, ket)` direction, the same thing
  `hexpr_pair_eval(drt, I, J, w)` returns. Wrap it with the dedup helper
  below, or call it in both orders, if you need G4's union.
- **It relies on G1.** The helper passes indices `0` and `1` into the
  sub-DRT it just built, taking row 0 to be the bra and row 1 the ket
  regardless of where the canonical walk places them. That is exactly
  the input-order guarantee, and it is why the helper is correct.

Per N4, the single-direction behavior of these helpers is current
behavior rather than contract. A future version may return the merged
decomposition, or add a separate entry point that does; code that wants
one specific direction today should call `pair_eval` with explicit
indices.

## Reconstructing the complete decomposition

To obtain the complete `(coef, pqrs)` decomposition of `H[I, J]`
from hexpr, call both directions and keep the larger of the two:

```python
def complete_pair_decomposition(drt, I, J, which="full"):
    """Return the complete per-pair (ij, pqrs, coef) decomposition
    of the matrix element H[I, J]."""
    a = drt.pair_eval(I, J, which=which)
    b = drt.pair_eval(J, I, which=which)
    return a if len(a[0]) >= len(b[0]) else b
```

(Equivalent helpers in Julia and Ruby follow the same pattern.)

A safer variant that does not rely on the "one direction contains
the other" property -- useful if one wants to defend against
future internal changes -- builds a `(ij, pqrs) -> coef`
dictionary from both directions:

```python
def complete_pair_dedup(drt, I, J, which="full"):
    ij_a, pqrs_a, coef_a = drt.pair_eval(I, J, which=which)
    ij_b, pqrs_b, coef_b = drt.pair_eval(J, I, which=which)
    merged = {}
    for ij, p, c in zip(ij_a, pqrs_a, coef_a):
        merged[(int(ij), int(p))] = float(c)
    for ij, p, c in zip(ij_b, pqrs_b, coef_b):
        merged.setdefault((int(ij), int(p)), float(c))
    # both forms work because coefficients agree on overlap
    return merged
```

The dedup form is recommended for new code: it does not depend on
the strict-subset property and adds an agreement check on shared
keys, so a future coefficient mismatch would surface loudly rather
than being silently absorbed. The keep-larger form is simpler and
correct on every case measured to date.

## Why hexpr keeps the directional split rather than merging it

It would be natural to ask why hexpr does not just emit the merged
decomposition itself. Two reasons.

First, the directional split is informative for studies of the
loop algebra. The "forward" and "backward" decompositions are the
two direction-selected traversals of the per-pair loop structure,
and a researcher studying that structure (or implementing a paper
such as Brooks's 1980 derivation) often wants to see each
separately. Merging them at the C level would discard that
distinction.

Second, the merge is cheap and local: a few lines of
binding-level code, easy to inline at the call site, easy to defer
until needed. The hexpr core stays a faithful realization of the
loop algebra; the user composes whatever view of the result they
need.

## API consequences for users

For CI Hamiltonian construction (the most common use case), users
who only need `H[I, J]` as a sum should either:

1. Wrap `pair_eval` as shown above and use the wrapper, or
2. Call `hexpr_subspace_eval` over the CSF set, which evaluates both
   orderings of every pair inside the set and is therefore complete.
   Note that this is a property of the subspace form specifically; a
   rectangular `hexpr_block_eval` with disjoint bra and ket sets is not
   complete. See "What block and subspace do, exactly" above.

For studies of the loop algebra itself (paper figures, debugging,
checking against alternative formulations), users can call
`pair_eval` directly in each direction and examine the two
decompositions side by side. This is the regime where the
directional split is informative.

## Relationship to the reference data (RubyRef)

The **complete** decomposition — the union of both directions —
matches the original reference data. The reference data in
`reference_pair/<case>/pair.txt` (F. Sasaki, T. Noro; the same
data the CI energies are anchored against to ~1e-8 Ha vs
OpenMolcas) stores both directions separately:

```
       I       J       pqrs     Coef
       1       2             4     -1.4142135623730949e+00
       2       1             4     -1.4142135623730945e+00
       2       1            53      1.4142135623730945e+00
       2       1            59     -2.8284271247461890e+00
       2       1            29     -1.4142135623730943e+00
```

`I < J` rows hold the forward decomposition; `I > J` rows hold the
backward decomposition. Cross-checking hexpr's per-pair union
against this file over all 190 pairs of s_full_4o4e gives 163
matches, 0 mismatches, 27 both-empty (163 = 137 + 26) — i.e. the
per-pair union reproduces the reference completely. The
per-direction *split* of that union (which terms land in the
forward vs backward set) is produced by hexpr's direction filter
as described above; the **physics** (the summed matrix element /
energy) is direction-independent and is the anchored quantity.

## Cross-references

- `c/include/hexpr.h` declares the C-level signatures.
- `c/README.md` links here from its "Note: internal sub-DRT API"
  section, which also documents the input-order guarantee G1 rests on.
- `docs/DECISIONS/directional-pair-decomposition.md` records the
  measurement that characterized this behavior and the decision
  not to alter the per-pair evaluator.
- `docs/PQRS_INDEX_SPEC.md` is normative for the `ij` and `pqrs`
  encodings referred to in G2.
- The Python / Julia / Ruby high-level bindings expose
  `csf_pair_eval(bra_steps, ket_steps)`; see "Step-based helpers in the
  bindings" above for what it returns and what it assumes.
