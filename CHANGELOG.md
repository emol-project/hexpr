# Changelog

All notable user-facing changes to hexpr are recorded here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/); versions follow
`HEXPR_VERSION_MAJOR.MINOR.PATCH` in `c/include/hexpr.h`.

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
