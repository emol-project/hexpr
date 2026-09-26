# Directional per-pair decomposition: mechanism, and why the per-pair evaluator was not rewritten

Status: **Accepted.**
Scope: the internal per-pair / block / subspace evaluator (`hexpr_pair_eval` and friends), its
directional decomposition behavior, and a proposal — not adopted — to replace its internals with a
DRT-free direct enumerator.
See also: `docs/PAIR_API_SEMANTICS.md` (user-facing semantics).

## Context

The per-pair evaluator, for an off-diagonal pair `(I, J)`, emits *one direction-selected*
decomposition of the symmetric matrix element `H[I,J]`, and `pair_eval(I,J)` vs `pair_eval(J,I)`
generally return different term sets (one a strict subset of the other). Internally each pair is
evaluated by building a small per-pair sub-DRT (`hexpr_drt_from_csfset`) and running the GUGA
loop walk on it.

Two questions were raised:

1. **Performance / honesty of the implementation.** Concern that advertising "evaluate a CSF *pair*"
   while internally *rebuilding a GUGA sub-DRT per pair* overstates the machinery — motivating a
   proposed rewrite to a "DRT-free direct per-orbital (`L_k`) enumerator."
2. **Is the directional asymmetry a feature or a defect?** The per-direction strict-subset structure
   looked like it might be an artifact of an internal deduplication step rather than intended
   behavior.

Both were investigated with read-only measurements before any code was changed.

## What was measured

**(a) `from_csfset` is already the minimal lattice.** The per-pair `hexpr_drt_from_csfset` only
inserts the `(orb,nve,is)` node codes visited by the two input walks (at most two per orbital level,
so `≤ 2(norb+1)` nodes), lex-sorts them, wires arcs, and runs the standard arc-weight / path-weight
passes. It does **not** run the heavy SOCI/FOCI construction (`hexpr_drt_create`: valence-node
generation, `prune_unreachable`, canonical-superset enumeration). The per-pair cost is `O(norb)` on a
tiny lattice, not a full DRT build.

**(b) The directional split is produced by the direction filter, not by deduplication.** In
`add_expr_pair` (`c/src/csf/csf_pair_core.c`), `delta = bra - ket` and loop elements with
`weight_g - weight_f != delta` are rejected. A controlled counterfactual on s_full_4o4e (190
off-diagonal pairs) attributed the behavior:

| variant | direction filter | loop-dedup (`rhash`) | fwd⊂bwd | bwd⊂fwd | both-zero | equal(>0) | incomparable | overlap-coef mismatch |
|---|---|---|---|---|---|---|---|---|
| baseline | on | on | 137 | 26 | 27 | 0 | 0 | 0 |
| rhash off | on | **off** | 137 | 25 | 27 | 1 | 0 | 0 |
| filter off | **off** | on | 0 | 0 | 0 | 1 | 189 | 603 |
| both off | **off** | **off** | 0 | 0 | 0 | 1 | 189 | 602 |

- The strict-subset shape (`incomparable = 0`, one side always ⊆ the other) is **present with `rhash`
  disabled** (137 unchanged). `rhash` is a minor loop-closure deduplicator; its only measured effect
  here is folding a single orientation pair. It is **not** the cause of the asymmetry.
- **Disabling the direction filter** does not symmetrize the output — it destroys the table (189/190
  pairs incomparable) and breaks coefficient agreement (up to 3.0). The filter is load-bearing, with
  direction also carried by the `bra`-dependent path-weight selection (`target_sum = bra - weight_g`)
  and the coefficient formula downstream of it.
- The "coefficients agree on overlapping keys" property follows from the direction filter (holds with
  `rhash` off; breaks only when the filter is removed).

**(c) The union is physically complete and matches the reference.** The per-pair union of both
directions matches `reference_pair/<case>/pair.txt` (Sasaki/Noro; the data the CI energies are
anchored against to ~1e-8 Ha vs OpenMolcas): over the 190 pairs, 163 match / 0 mismatch / 27
both-empty. The summed matrix element (the physics) is direction-independent.

**(d) The directional contract is intrinsically a sub-DRT property.** A DRT-free enumerator that
walks only the two fixed paths reproduces the per-orbital `∏ L_k` values bit-exactly (validated as an
independent cross-check), but it emits the *symmetric* full set in both directions. Reproducing the
current per-direction contract byte-exactly requires the sub-DRT's integer weight bookkeeping
(`weight_g`/`weight_f` and the orientation handling) — i.e. exactly what `from_csfset` already
computes. A parallel old-vs-new harness, run over every pair (not only the frozen regression rows),
surfaced this directly on the first spin-recoupling pair.

## Decision

**Do not rewrite the internal per-pair evaluator.** `pair_via_subdrt + from_csfset` is already the
minimal-lattice-plus-kernel design the rewrite aimed for; `from_csfset` is the `O(norb)` two-walk
lattice, and the directional contract is defined by the integer weight bookkeeping it produces.
A byte-exact "DRT-free" rewrite would have to reconstruct that same bookkeeping for **no correctness
gain and nonzero risk**; the only thing it could shed is a small per-pair allocation, which is a
performance concern explicitly out of scope.

The directional asymmetry is retained as the documented behavior. It is a property of the direction
filter / GUGA loop walk, not a deduplication artifact; the physics (summed element) is
direction-independent and anchored; and consumers already obtain the complete decomposition by
merging both directions (`docs/PAIR_API_SEMANTICS.md`).

## Consequences

- **No change to public API, IO, index conventions, or `walk_*`.** The frozen regression oracle
  (`reference_csfset/`, `c/tests/test_csfset_pair.c`) passes unchanged (37/37). After reverting the
  exploratory scaffolding, the library's exported symbol set is identical to the frozen baseline.
- The DRT-free value-path enumerator was kept as an **independent `∏ L_k` cross-check**, out of the
  library build (`c/scratch/`); it is committed for reference but is **not** part of the shipped
  library and is stripped from the distribution before public release (with the other pre-release
  scaffolding).
- **Open, by design:** *why* the surviving structure is always a strict subset rather than merely
  overlapping is an analytic property of the GUGA loop algebra (reduction / telescope; see the
  methodology paper), not established by the build-level counterfactual above, which only shows the
  property is direction-filter-governed and `rhash`-independent.
- **If revisited:** emitting the symmetric full decomposition per direction (instead of the
  direction-selected split) is a *contract change* — it requires re-freezing the per-pair oracle and
  assessing the one in-tree consumer that depends on the merged-output contract — and would need
  explicit sign-off. Not chosen, because the physics is already correct and direction-independent and
  downstream consumers merge both directions.

## Provenance note (verify before publish)

The frozen regression artifacts correspond to the pre-release tree, not to the current branch HEAD.
The relationship, confirmed from the commit history:

- **`cb8a3f3`** — the commit at which the CSFSET oracle *data* was frozen, per
  `reference_csfset/README` ("the EXACT output the CURRENT (commit cb8a3f3) sub-DRT path produces";
  it established the input-order directional contract). This is the regression baseline.
- **`bc2c6fd`** — the tree from which the frozen `libhexpr.so` binary (mtime 2026-06-18 23:14) was
  built; it follows `cb8a3f3`.
- **`b492548`** — current HEAD prior to the docs work, adding only the warning-hardening (Topic 1)
  and macOS-support (Topic 2) commits on top — none of which touch the directional core.

Ordering is `cb8a3f3 → bc2c6fd → b492548`, with no change to the directional core
(`from_csfset`/arc-weights, the loop walk and `rhash`, the direction filter, the 16 operator cases,
`add_expr_pair`) at any point — verified byte-identical across these commits (Topic-3 integrity recon
+ the directional-cause probe).
