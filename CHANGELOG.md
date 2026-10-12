# Changelog

All notable user-facing changes to hexpr are recorded here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/); versions follow
`HEXPR_VERSION_MAJOR.MINOR.PATCH` in `c/include/hexpr.h`.

## [1.1.0] -- 2026-10-12

### Added

- Complete-decomposition calls: `hexpr_pair_count_full`,
  `hexpr_pair_eval_full`, `hexpr_block_count_full` and
  `hexpr_block_eval_full`, with their F77 wrappers (`hexpr_pair_count_full_`
  and so on) and the corresponding functions in the Fortran module and the
  Python, Julia and Ruby bindings (`pair_count_full`, `pair_eval_full`,
  `block_count_full`, `block_eval_full`, `csf_pair_eval_full`).
  `hexpr_pair_eval` returns one direction-selected decomposition of
  `H[I,J]`, and the complete one is the union of both argument orders
  (`docs/PAIR_API_SEMANTICS.md`, G4). The new calls return that union in one
  call, the same terms with the same coefficients for either argument order,
  and each `(ij, pqrs)` key once. `hexpr_block_eval_full` is therefore
  complete for a block whose bra and ket sets are disjoint, which
  `hexpr_block_eval` is not.
- They are also faster than calling `hexpr_pair_eval` in both orders: the
  loop search follows only the bra and ket walks. Measured per pair on
  random pairs, against two `hexpr_pair_eval` calls: 56 vs 216 us (CAS(6,6)
  singlet, 175 CSFs), 59 vs 207 us (FOCI, 125 CSFs), 65 vs 345 us (6
  electrons in 8 orbitals, 1176 CSFs).
- `c/tests/test_pair_full_union.c` checks, on every pair of SOCI, FOCI,
  full-CAS (singlet, triplet, quintet) and CSFSET spaces, that the new calls
  equal the union of the two directions of `hexpr_pair_eval` with agreeing
  coefficients, independent of the argument order, for `HEXPR_EXPR_FULL` and
  `HEXPR_EXPR_ONE`, and that the block form matches the pair form.

### Fixed

- Pair and block evaluations can now run concurrently in several threads,
  on one shared DRT handle or on separate ones. In 1.0.0 they shared a
  process-wide L(k) cache (filled on a miss), and the README's statement
  that read-only queries on a handle are reentrant did not hold for them.
  The cache and the walk restriction of the `*_full` calls are now
  thread-local (C11 `_Thread_local`), and the sqrt(2), sqrt(3) table is a
  constant. `c/tests/test_pair_threads.c` checks that 8 threads give
  bit-identical terms to a serial run for `hexpr_pair_eval` and
  `hexpr_pair_eval_full`, on a shared handle and on per-thread CSFSET
  sub-DRTs. The serial speed is unchanged (the search fetches the
  thread-local state once per call and passes it down).

### Changed

- Faster pair evaluation, same output. The loop-deduplication table of the
  search was cleared (1,024 buckets) at every start node and walked in
  full when released; it now records the buckets it used and resets only
  those. The construction tables of a CSFSET DRT are sized to the set
  instead of 4,096 buckets, and the search keeps its per-level arrays on
  the stack up to 128 orbitals. On the pairs of a real MRCI build (N2,
  10 electrons in 8 internal orbitals) one complete pair decomposition
  takes 57 instead of 93 us including the sub-DRT; with many threads the
  gain is larger (less allocator and memory traffic), e.g. 36.6 -> 12.7 s
  for a 37,162-CSF compressed build on 16 threads in hexpr2. Every
  test, including the byte-identical frozen oracles, passes unchanged.

### Unchanged

- Every v1.0 call returns exactly what it returned in 1.0.0; the frozen
  `reference_csfset/` oracle is reproduced byte-for-byte. The new calls do
  not alter the directional contract of the existing ones.

## [1.0.0] -- 2026-09-26

Initial public release.

### Known limitations

- `hexpr_drt_create` does not validate `nve` against the
  `(ncore,nval,next)` partition in general; of the inconsistent
  combinations it rejects only `ncore=0` together with `nval=0`. See the
  "still open" section in `docs/DRT_CONSISTENCY_INVESTIGATION.md`.
- On macOS, `libhexpr.dylib` also exports the five `wigner_*` entry
  points of the bundled Wigner library. They are hidden on Linux with
  `--exclude-libs`, which has no Darwin equivalent. They are not part of
  the API and may disappear; link `libwigner.a` directly if you need
  them (see `c/wig/README.md`).

### Notes on pre-publication changes

Before this release the code was used internally as an unpublished
build. Changes made since then, recorded here because they define
current behaviour:

- `libhexpr.so` exports only the API: the 45 functions declared in
  `c/include/hexpr.h`, the 45 F77 wrappers used by the Fortran binding,
  and `hexpr_fortran_stdout`. It previously exported 142 symbols -- the
  rest being internal entry points that no installed header declares,
  among them generically named ones such as `darr_push` and `dhash_get`,
  which an application defining the same name could interpose. Code that
  reached an undeclared symbol will no longer link; there was never a
  declaration for it to use. The Wigner routines remain available on
  their own terms through `libwigner.a` (`c/wig/README.md`). No
  declaration in `hexpr.h` changed, and the ABI of every declared
  function is unchanged.
- `make -C c` now also builds `c/libhexpr.a`, a self-contained static
  archive over the same objects. Programs that legitimately need an
  internal entry point -- the in-tree test suite and `run_reference` --
  link that instead of the shared library.
- `hexpr_pair_count/_eval`, `hexpr_block_count/_eval` and
  `hexpr_subspace_count/_eval` take INPUT-ORDER CSF indices, in
  `[0, hexpr_drt_ncsfs())`, the same space as `hexpr_drt_csfs`,
  `hexpr_csf_index` and `hexpr_csf_steps`. They previously took canonical
  (walk-order) indices, which differ from input order for a CSFSET
  sub-DRT: a loop over `0..hexpr_drt_ncsfs()-1` was writable, correct for
  SOCI/FOCI, and silently returned matrix elements of other CSFs for
  CSFSET, since the indices stayed in range. There is now one index space
  across the C API and all five bindings. The emitted `ij` is unchanged
  and remains a canonical pair address per `docs/PQRS_INDEX_SPEC.md`; use
  `bra_pos`/`ket_pos` to attribute a term to a pair. No declaration in
  `hexpr.h` changed, and the frozen `reference_csfset/` oracle is
  reproduced byte-for-byte.
- SOCI/FOCI DRTs: `hexpr_drt_ncsfs()` / `hexpr_drt_walk_ncsfs()` now always
  agree with the number of CSFs actually enumerated (an earlier over-count
  could cause out-of-bounds reads and silently empty pair results). See
  `docs/DRT_CONSISTENCY_INVESTIGATION.md`.
- `hexpr_drt_create` and `hexpr_drt_from_csfset` return `NULL` (with
  `hexpr_set_error`) if the enumerated CSF count disagrees with the
  reported count, or if the first node is not the walk origin, instead of
  returning an inconsistent DRT.
- The expression-generation kernel is table-driven (18-entry one-orbital
  matrix-element table, a declarative Brooks-case table shared by the
  whole-matrix and per-pair drivers, precomputed per-orbital factors L(k)).
  Generated expressions are byte-identical on the full reference suite.
