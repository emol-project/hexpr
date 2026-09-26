# DRT consistency investigation: CSF-count divergence, dead nodes, and what it did and did not turn out to be

> Note: commit hashes in this document refer to the private development
> history and do not exist in the public repository, which starts from a
> single snapshot commit.

Status: **Resolved** for the two root causes found and fixed (see below); **partially open** on
three named points (nve-vs-partition validation, the `hexpr_drt_ncsfs`/`hexpr_drt_walk_ncsfs`
naming duplication, and a top-is-last sort assumption). This document is the full account, written
for someone with no other context. It records the false starts as they happened, not only the
conclusions that survived them.

See also: `CHANGELOG.md` (user-facing summary), `c/README.md` (the corrected CSF-ordering
guarantee), `docs/PAIR_API_SEMANTICS.md` (a related, separate, and correct piece of behavior that
came up during this investigation and was cleared of suspicion).

## How this started

This did not begin as a hexpr investigation. It began downstream, in **hexpr2**, a separate
project that uses hexpr as its evaluation engine. hexpr2's `hexpr2.CsfSet.classified` constructor
calls `hexpr_drt_create` directly (it is the only one of hexpr2's four CSF-set constructors that
does -- `distribution`, `internaln`, and `excitation` build CSFs via their own pure-integer code
and never touch `hexpr_drt_create`). While choosing parameters for an unrelated invariant sweep in
hexpr2's own test matrix, a specific `(kind, nve, spin_x2, ncore, nval, next)` combination was
checked for exposure to a shape of defect already on record in that project's notes, and it landed
on a reproducible problem:

```
>>> hexpr2.CsfSet.classified(hexpr2.SOCI, nve=6, spin_x2=2, ncore=1, nval=1, next=2)
ncsfs= 2
[[  3   3   1   1]
 [103 127   0   0]]     # run 1
[[  3   3   1   1]
 [ 68 127   0   0]]     # run 2
[[  3   3   1   1]
 [221 127   0   0]]     # run 3
```

Only one real CSF (`[3,3,1,1]`) exists; `ncsfs` reported 2, and the second row was not zero or a
fixed wrong value -- it varied between separate process runs, with one byte position stable and the
other not. That is the signature of reading memory the allocator gave back but the library never
wrote into, not a fixed logic error. This is why the first hypothesis, below, was memory
corruption: that is what the evidence looked like on first contact, from outside hexpr entirely.

Reproducing the same call directly through hexpr's own Python binding (bypassing hexpr2 completely)
showed the identical symptom:

```python
d = hexpr.Drt(hexpr.HEXPR_DRT_SOCI, nve=6, spin_x2=2, ncore=1, nval=1, next_=2)
d.ncsfs          # 2
d.csfs()
# [[  3   3   1   1]
#  [131 127   0   0]]      # garbage, varies across runs, same 127 tail byte
```

which located the problem inside hexpr itself (`hexpr_drt_create`/`hexpr_drt_ncsfs`/
`hexpr_drt_csfs`), not in hexpr2's own code, and moved the investigation into this repository.

A reader who assumes this was found by fuzzing hexpr, or by hexpr's own test suite, would misjudge
what the actual coverage gap was: it is not that hexpr lacked tests: at the time, no test in this
repository, and no shipped example, happened to choose parameters landing on the trigger (see
"What was verified unaffected" below for the specific check of that claim). It surfaced because a
downstream consumer's test matrix, built for its own reasons, happened to probe the right corner.

## The chain of hypotheses, in the order they were raised

Five explanations were considered over the course of this investigation. Three held up; two were
raised, investigated, and then explicitly retracted.

1. **Memory corruption.** *(Superseded, not wrong exactly -- refined.)* The garbage bytes'
   run-to-run variation is real and is explained, but not by corruption in the sense of something
   overwriting live data: it's an **out-of-bounds read**. `hexpr_drt_csfs()`'s buffer is allocated
   to hold exactly as many rows as `hexpr_drt_mk_csf_list` actually enumerated (3, in the example
   above); `hexpr_drt_ncsfs()`/`hexpr_drt_walk_ncsfs()` reported 5 for the same DRT in the
   regression case ultimately adopted (`SOCI, nve=6, spin_x2=0, ncore=1, nval=1, next=2` -- a
   smaller, cleaner reproduction than the `spin_x2=2` case above, found during the investigation
   and used for everything from here on). A caller sizing a read by the reported count reads past
   the allocation. Confirmed directly with Valgrind:
   ```
   ==4== Invalid read of size 1
   ==4==    at 0x1092E1: main
   ==4==  Address 0x4b984bc is 0 bytes after a block of size 12 alloc'd
   ==4==    at 0x4844818: malloc (vg_replace_malloc.c:446)
   ==4==    by 0x4860F4F: hexpr_drt_create (drt.c:1223)
   ```
   12 bytes = 3 rows x 4 orbitals, the true count; the reported count said 5. Every language
   binding's `.csfs()`/`.walk_csfs()` sizes its read this same way (checked directly: Fortran,
   Julia, Ruby, Python all size by `ncsfs()`/`walk_ncsfs()`), so this was not a Python-specific or
   hexpr2-specific issue.

2. **An allocation bound reported as an exact count.** *(Also superseded, by finding it isn't even
   reliably that.)* Once the OOB read was understood, the next candidate was that
   `hexpr_drt_ncsfs()`/`hexpr_drt_walk_ncsfs()` were meant as an upper bound on allocation size, not
   an exact row count, and callers were simply misusing them. This was checked against both the
   header documentation and the implementation: `hexpr.h` documents both functions as exact counts
   (for SOCI/FOCI, explicitly "equals `hexpr_drt_csfs()`/`hexpr_drt_ncsfs()`", never described as a
   bound), and the implementation backs that reading -- both functions return the identical field
   (`nodes[0].number_upper_arcs`), computed by `store_number_upper_arcs`'s recursive
   `num += imax(child_n, 1)`, a floor added only to credit the top node's own single walk. In a
   graph where every non-top node has a genuine path to the top, this floor is a no-op and the
   count is exact. It is not, however, derived by counting completed walks -- so it is only
   *conditionally* exact, and nothing documents that condition or checks it. This ruled out reading
   the field as an intentional bound: it is presented, and normally behaves, as an exact count, and
   the actual defect is that the condition for exactness can silently fail.

3. **The walk enumerator generating invalid walks.** *(Investigated, not confirmed -- the opposite
   direction was checked instead, and found clean within what was tested.)* The concrete concern
   was the mirror-image risk: given that the count can be *too high*, could the DFS walk enumerator
   (`hexpr_drt_mk_csf_list`/`tree_search_csf`) itself ever be *too low* -- silently dropping CSFs
   that a correct GUGA/Paldus admissible space should contain? An over-count is visible (empty rows
   in a Hamiltonian); an under-count would not be, since the shipped examples only take the lowest
   eigenvalue. This was checked, not proven in general: on eight `(n, N, S)` configurations with
   `ncore=0, next=0` (chosen because the Paldus/Weyl dimension formula gives an independent exact
   count there), `hexpr_drt_ncsfs()` matched the formula exactly in all eight, and separately, the
   DFS itself was shown to be structurally incapable of *missing* an existing arc (it is an
   unconditional traversal of every `upper_arc` actually present, with no pruning of its own) --
   only `store_number_upper_arcs`'s floor can inflate, never deflate. This is evidence that the
   defect's direction is over-counting, not proof that no under-counting configuration exists
   anywhere; it was not chased further because nothing in the rest of the investigation ever
   pointed the other way.

4. **Wrong-`ij` attribution in `hexpr_expr_build`.** *(Retracted -- not a hexpr defect, but a real
   question worth asking, and the retraction is clean.)* If the phantom node's inflated
   `number_upper_arcs` shifts `arc_weight` prefix sums for its siblings, could a *real* term end up
   labeled with the wrong canonical pair index (`ij`), rather than merely being missing? This was
   tested directly, not inferred: `hexpr_expr_build`'s full-DRT output was diffed term-by-term
   against a pairwise ground truth built from isolated `hexpr_drt_from_csfset` pairs, on both the
   pathological DRT and a healthy control. Result: **zero terms under the wrong `ij`**, in both
   cases -- every term `hexpr_expr_build` emitted matched an entry in the ground truth at the same
   `ij`. This hypothesis is cleanly false, and the investigation said so once it had the diff to
   show, rather than continuing to suspect it.

5. **An off-diagonal term deficit.** *(Retracted in full, and the retraction is the important part
   of this whole document.)* Diffing `hexpr_expr_build` against the pairwise ground truth (for
   hypothesis 4, immediately above) surfaced a real-looking count discrepancy: 55 terms from the
   full build vs. 139 from pairwise reconstruction on the pathological DRT, with all 84 missing
   entries concentrated on the off-diagonal pairs. A control on a completely healthy 3-CSF DRT (no
   phantom node possible) reproduced the identical shape -- same kind of deficit, same concentration
   on off-diagonal pairs -- which meant it could not be caused by the phantom-node bug at all. At
   that point there were two candidate explanations still open: (a) the deficit is a real defect in
   `hexpr_expr_build`, or (b) it is an artifact of how the pairwise "ground truth" was assembled. Both
   were pursued, and (b) is what turned out to be true; see the next section for how that was
   settled and what the mistake actually was.

**Both retractions (4 and 5) were self-corrections, not conclusions handed down from re-reading
hexpr's source more carefully.** Hypothesis 4 was settled by an actual diff of emitted terms, not by
re-inspecting the `arc_weight` code path a second time. Hypothesis 5's retraction came from finding
a bug in the comparison harness itself, confirmed independently (below) -- not from hexpr's code
being re-read and cleared. A reader who sees only "wrong-ij: not confirmed" and "off-diagonal
deficit: not a real defect" without this paragraph would have no way to tell how much of this
account is hexpr's code holding up under test versus the investigation correcting its own tools;
that distinction matters for how much confidence to place in the parts that did stand (case 4
survived a direct term-by-term diff; case 5's surviving conclusion is corroborated by an
independent physical calculation, described next, not merely by re-reading a fixed comparison
script).

## What made each hypothesis decidable

Hypotheses 1-3 were settled by direct instrumentation: Valgrind for the OOB read (1), reading the
header contract against the implementation for the bound question (2), and a formula-based ground
truth (Paldus/Weyl dimension) plus structural argument about the DFS for the undercount question
(3).

Hypothesis 5 (the off-diagonal deficit) is the one that most needed an authority independent of
hexpr's own code, because by that point the investigation's own comparison harness was itself a
suspect, and re-reading it is exactly the kind of self-check that can miss a bug in its own
premises. It was settled by building the full Ms=0 Hamiltonian for the same arbitrary
`h1`/two-electron-integral inputs directly through PySCF's `fci.direct_spin1.contract_2e` -- mature,
independently-validated code with no dependency on hexpr at all -- and diagonalizing it:

```
FCI (independent, 4 determinants): [-0.594262, -0.300000, 1.642347, 4.451915]
H_A (hexpr, whole-DRT build):      [-0.594262,             1.642347, 4.451915]   <- 3 of the 4, exact
H_B ("ground truth" harness):      [-4.04,      0.32,      9.22               ]  <- matches nothing
```

`H_A` -- built the whole-DRT way, the thing under suspicion -- reproduced **three of the four**
independent FCI eigenvalues to six decimal places. The fourth FCI eigenvalue, `-0.300000`, is
exactly the Ms=0 **triplet** state that a 3-CSF **singlet** DRT has no business containing: the FCI
calculation works in the full 4-determinant Ms=0 space, which includes one triplet component
alongside the three singlet states; hexpr's own DRT is spin-adapted and singlet-only, so its
absence from `H_A`'s spectrum is correct, not a deficiency. `H_B`, the harness's own "pairwise
ground truth," matched none of the four independent eigenvalues -- it was `H_B`, not hexpr, that was
wrong.

The actual harness bug: reading `hexpr_expr_build`'s output off an isolated 2-CSF
`from_csfset(a,b)` sub-DRT and reading its *entire* term list (each CSF's own diagonal plus the
genuine cross terms) as if it were purely the off-diagonal coupling, without decoding each term's
local `ij` and discarding the local self-terms first. That double-counted each CSF's diagonal
energy into the off-diagonal bucket, inflating it and producing the 84-term "deficit" against the
correctly-smaller true cross-coupling. Once corrected, `H_A == H_B` to machine precision on the
healthy control, and re-running the original 55-vs-139 comparison on the pathological DRT gave 55
of 55 exact, zero missing, zero extra. `hexpr_expr_build` was never implicated; this is now noted
as a caveat for `docs/PAIR_API_SEMANTICS.md`'s existing recommendation to use
`from_csfset(pair).terms_full()`: it returns the *whole* expression for that sub-DRT, not only the
off-diagonal coupling, and a caller reconstructing a multi-CSF matrix from it must decode and filter
accordingly.

## Root cause, the two fixes, and the assertions added

The root cause, in both its forms, lives in `prune_unreachable` (`c/src/loopgen/drt.c`), the
bottom-to-top / top-to-bottom graph reachability sweep that decides which generated nodes survive
into the final DRT. `generate_valence_node` inserts every node it builds into the construction
table regardless of its own local admissibility (`status` is set independently, but an invalid node
is still inserted, not discarded). `prune_unreachable`'s two BFS passes then need to treat an
invalid node as radiating no reachability, in **both** of the ways that could matter, and the code
only did one of them:

1. **A neighbor's own invalidity was not checked before crossing into it** (fixed in `a6acbf4`).
   A node reachable only through an already-invalid intermediate was marked reachable anyway, then
   survived into the final sorted `nodes[]` array with no real arc back through that (excluded)
   intermediate once arc-wiring ran (arcs are only wired among survivors). That "dead node" has an
   incoming arc but no real outgoing path to the top; `store_number_upper_arcs`'s
   `imax(child_n, 1)` floor then credited its parent with one phantom walk, inflating
   `nodes[0].number_upper_arcs` -- and hence both `hexpr_drt_ncsfs()` and `hexpr_drt_walk_ncsfs()`
   -- past what `hexpr_drt_mk_csf_list` actually enumerated and past what `csf_packed` was
   allocated to hold. This is the cause behind hypotheses 1 and 2 above, and behind the original
   `nve=6, ncore=1, nval=1, next=2` reproduction.

2. **The dequeued node's own invalidity was not checked before expanding its successors/
   predecessors** (fixed later, in `25e1075`, in a separate round of work on this same
   investigation). This is the same omission, symmetrically: a node's own `status==0` should mean
   it radiates no reachability of its own, not merely that it fails to be crossed *into*. It went
   unnoticed at first because it only matters when the very *seed* of a BFS pass -- the bottom
   `(0,0,0)` or the requested top state -- is itself the invalid node. That happens specifically
   when a construction has `ncore=0` and `nval=0`: `generate_valence_node`'s own admissibility bound
   makes the origin trivially inadmissible whenever there is no core or valence reference to supply
   any of the requested electrons and the requested `nve` exceeds SOCI's <=2-particle / FOCI's
   <=1-particle external-region budget. Before this second fix, the forward BFS still seeded
   reachability at that invalid bottom and expanded through it anyway, so a spurious path from a
   neighboring node could still reach the top, and `hexpr_drt_create` would return a DRT built one
   orbital short of the true origin instead of correctly finding that no valid walk exists at all
   and rejecting the request. **This is not a special case for `nval=0`** -- the fix
   (`if (!cur->status) continue;` in both BFS passes, symmetric with the existing neighbor check)
   is the same general graph-reachability rule applied to the node actually being expanded, not
   only to the nodes it might expand into. `nval=0` is simply the smallest combination that
   exercises it; the general predicate `hexpr_drt_create` enforces (unaffected by either fix) is
   "there exists at least one bottom-to-top walk consistent with the kind's step rules and
   region budgets," and construction should succeed whenever that holds and fail when it does not,
   regardless of which region is empty.

Alongside the two graph fixes:

- **A construction-time consistency guard** (`7101e11`): both `hexpr_drt_create` and
  `hexpr_drt_from_csfset` now compare the `mk_csf_list`-enumerated count against
  `nodes[0].number_upper_arcs` immediately after building the CSF list, and reject construction
  (`hexpr_set_error` + `NULL`) on any mismatch, rather than publishing an inconsistent DRT. This
  guard is blind to the `nval=0` cause on its own: `mk_csf_list` and `store_number_upper_arcs` both
  start from the *same* `nodes[0]`, so if it is the wrong node, both sides shift together and still
  agree with each other -- which is exactly why the second fix and the invariant below were still
  needed after this guard existed.
- **An explicit invariant that `nodes[0]` is the walk origin `(0,0,0)`** (`c5ba528`), placed next to
  the guard above, in both constructors. Checked directly, not assumed: with the neighbor-status fix
  (`a6acbf4`) applied but the self-status fix (`25e1075`) reverted, this assertion alone catches all
  four `nval=0` trigger combinations (`nodes[0]` comes back as `(orb=1,...)`, not `(0,0,0)`). With
  both fixes in place, it does not fire on any of the 970 swept combinations -- the pre-existing
  "no node exists for the requested target state" check now catches those four cases earlier in the
  pipeline. It remains as defense in depth for any future or unforeseen way `nodes[0]` could end up
  wrong, the same role the count-mismatch guard already plays for the count divergence.
- **Two status discards closed in `csf_pair.c`** (`066498f`): `pair_via_subdrt` was calling
  `hexpr_drt_walk_csf_steps`/`hexpr_drt_walk_csf_index` and discarding the returned status. Under
  Valgrind, a decode failure on the pathological DRT was shown to write only part of the output
  buffer and leave the rest as uninitialized memory, which was then fed into
  `hexpr_drt_from_csfset` -- silently returning zero terms if the uninitialized byte happened to
  fall out of `0..3`, or, worse, silently building a sub-DRT for a fabricated wrong CSF if it
  happened to fall in range. Both call sites now check the status and return empty on failure.

Commit sequence (`release/v1.0.0`, oldest first): `b1013b8` (regression test for the first cause),
`7101e11` (count-mismatch guard), `a6acbf4` (neighbor-status fix, the first root cause),
`066498f` (status-discard cleanup), `6b980f5` (README guarantee), `764b64e` (regression test for
the second cause), `25e1075` (self-status fix, the second root cause), `c5ba528` (nodes[0]
invariant), `0a604dc` (top-is-last note, below).

## What is still open

Three things are named here deliberately, not fixed, and not to be read as resolved by anything
above:

1. **`hexpr_drt_create` does not validate `nve` against the `(ncore, nval, next)` partition in
   general.** It rejects a request only if no node exists for the resulting target state after
   generation and pruning (and, as of the `ncore=nval=0` fix, correctly does so for that specific
   shape). There is no check relating `nve`/`spin_x2` to region capacity up front, and no comment
   anywhere asserting an intended range. Attempts during this investigation to derive a general
   `ncore`-exclusion capacity formula from `hexpr.h` did not succeed (varying `ncore` while holding
   `nval+next` fixed gave a non-monotonic accept/reject pattern inconsistent with a simple
   "`nve` excludes core" bound) -- consistent with SOCI/FOCI's admission rule being the RAS-style
   "<=2 holes / <=2 particles" interacting-space definition mentioned elsewhere in this repository's
   docs, not a flat capacity bound. That formula is not written down anywhere in code or docs; this
   investigation did not reverse-engineer it.
2. **`hexpr_drt_ncsfs()` and `hexpr_drt_walk_ncsfs()` are two names for one computation, for
   SOCI/FOCI.** Both return the identical field (`nodes[0].number_upper_arcs`) for every DRT kind
   except `CSFSET`, where they differ (input count vs. canonical superset count). `hexpr.h`
   documents this equality as intentional, but it means the two-name distinction only carries
   meaning for one of the three DRT kinds; whether that is the right shape for the public API, or
   whether SOCI/FOCI should have just one name, was raised during the fix-options discussion and
   not decided.
3. **The top-is-last sort assumption** (`0a604dc`): `hexpr_drt_create`'s
   `store_number_lower_arcs(nodes, valid_count - 1)` assumes the top node sorts last in the final
   `nodes[]` array (by `node_cmp`'s ascending `(orb,nve,is)` order) -- unchecked, the mirror image of
   the `nodes[0]`-is-bottom assumption this investigation did assert. `hexpr_drt_from_csfset`'s
   analogous call passes `top_index` explicitly instead of relying on sort position, so the two
   constructors are inconsistent with each other here. Recorded in the source where it lives; not
   fixed, since nothing in this investigation's testing (970-combination sweep included) found it to
   actually diverge from `top_index`.

## What was verified unaffected, and how

- **The three shipped `example_ci/h2o_soci` configurations** (`run_pyscf.py`, `run_psi4.py`: SOCI;
  `run_pyscf_foci.py`: FOCI; all `nve=8, spin_x2=0, ncore=2, nval=4, next=11`, `norb=17`). Checked
  directly against both the original (pre-investigation, `6ee64c8`) and the fully-fixed library:
  `ncsfs`, `walk_ncsfs`, and the `mk_csf_list`-enumerated count agree exactly in both configurations
  (SOCI: 17685; FOCI: 2140), a full `hexpr_drt_csfs()` buffer read is Valgrind-clean in both, and
  `hexpr_pair_eval(drt, 0, 0, ...)` matches `hexpr_drt_from_csfset`-based ground truth exactly for
  both (SOCI: 20 terms, coefficient sum 24; FOCI: 38 terms, coefficient sum 26). The full assembled
  calculation was then run end-to-end, original library vs. fixed library, and the results are
  byte-for-byte identical:

  | | ncsfs | nterms | Energy (Hartree) |
  |---|---|---|---|
  | SOCI, original library | 17685 | 21,437,024 | -76.1971171060 |
  | SOCI, fixed library | 17685 | 21,437,024 | -76.1971171060 |
  | FOCI, original library | 2140 | 1,256,191 | -76.0630642167 |
  | FOCI, fixed library | 2140 | 1,256,191 | -76.0630642167 |

  (Nuclear repulsion 9.1892096261 Ha, RHF -76.0091065449 Ha, H2O/6-31G(d), spherical basis, both
  runs identical there too.) This is confirmation, not a discovery: `nterms` for both kinds already
  agreed before this number was pasted here, and neither `hexpr_expr_build`/`brooks_all` nor the
  assembly scripts were touched by anything in this investigation -- the energy match is expected,
  not a surprise, and is recorded because a number cited in a paper should be seen run, not only
  inferred from unrelated evidence.

- **The downstream `hexpr2` case that started this investigation.** hexpr2's own `classified` case
  7 (`SOCI, nve=4, ncore=1, nval=1, next=2`, `xval/hexpr2_ch2_classified_1c1v2x_4e.py`) was
  re-checked directly, separately from the `nve=6`/`nve=6,spin_x2=2` cases that did trigger the
  divergence: 15 CSFs, 0 corrupted rows, byte-identical (SHA-256) across three separate process
  runs, both before this investigation's fixes existed and after. It was never in the trigger
  region (clean `nve<=4` for that `ncore/nval/next` slice); nothing in this investigation required
  or produced any change to hexpr2 itself.

- **Every existing test in this repository.** `c/tests` (`make check`), the Fortran suite, and the
  Python, Ruby, and Julia binding test suites all pass unchanged against the fixed library, with one
  addition: `c/tests/test_ncsfs_consistency.c`, the regression test this investigation added, now
  passes where it previously (and intentionally, for its first two cases) failed.
