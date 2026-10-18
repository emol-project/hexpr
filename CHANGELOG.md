# Changelog

All notable user-facing changes to hexpr are recorded here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/); versions follow
`HEXPR_VERSION_MAJOR.MINOR.PATCH` in `c/include/hexpr.h`.

## [1.2.0] -- 2026-10-18

### Added

- Spin-dependent one-electron operators: `hexpr_rank1_pair` and
  `hexpr_rank1_max_terms` return the reduced matrix elements
  `<S'||(a+_p a~_q)^1||S>` (Edmonds convention) between two CSFs given as
  step vectors, which may differ in total spin; no DRT is built. The
  per-orbital factors are those of the spin-free path; the operator rank is
  carried to the top orbital and the result is multiplied by
  `(-1)^(1+S'-S) sqrt(2S+1)`. F77 wrappers `hexpr_rank1_pair_` and
  `hexpr_rank1_max_terms_`, Python `hexpr.rank1_pair`, Fortran module
  `hexpr_rank1_pair` / `hexpr_rank1_max_terms`, Julia `HExpr.rank1_pair`,
  Ruby `HExpr.rank1_pair`, each with a test of the worked examples, the
  total-spin identity and the adjoint relation.
- `c/tests/test_rank1.c`: total-spin identity, adjoint relation over every
  CSF pair up to 5 orbitals, worked examples, selection rules. During
  development the formulation was checked against determinants for ranks
  0, 1, 2, and the C evaluator against determinants (max 4.4e-15).
- `example_ci/ch2_so`, `example_ci/oh_so`: spin-orbit matrices from PySCF
  integrals, checked against determinants and against GAMESS HSO1, HSO2P
  and HSO2 (the GAMESS inputs are included; its outputs are not).
- `validation/safety/`: compiler and sanitizer matrix (gcc/clang at
  -O0/-O2/-O3, ASan+UBSan, integer sanitizers, MSan, TSan), valgrind,
  static analysis (cppcheck, gcc `-fanalyzer`) and fuzzing of
  `hexpr_rank1_pair`, with the results for this release.
- `c/tests/test_oom.c` (Linux): fails each allocation of a call in turn and
  checks that the call reports out of memory or returns its normal result,
  with no state left behind; under ASan/LSan also no leak or invalid access.

### Fixed

- Out of memory is now reported instead of crashing or returning a wrong
  result. Before, a failed `malloc`, `calloc` or `realloc` inside DRT
  construction, expression generation or pair evaluation was dereferenced
  (a crash), and in the pair/block/subspace calls some failures returned
  `HEXPR_OK` with no terms; the stable sort fell back to an unstable `qsort`,
  which could reorder the output. Now the calls that return a status
  (`hexpr_pair_*`, `hexpr_block_*`, `hexpr_subspace_*`, `hexpr_csf_to_dets`)
  return `HEXPR_ERR_OOM`, and the builders (`hexpr_drt_create`,
  `hexpr_drt_from_csfset`, `hexpr_expr_build`) return NULL with "out of
  memory" in `hexpr_last_error()`; everything allocated so far is freed.
  A pair whose walk cannot be decoded (ruled out by the index validation)
  is reported as `HEXPR_ERR_INTERNAL` instead of as no terms.

### Unchanged

- Every 1.0.0 and 1.1.0 call and its output when memory suffices (the
  frozen references in `reference/`, `reference_pair/`, `reference_csfset/`
  pass unchanged).

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
