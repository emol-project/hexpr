# validation/safety -- memory-safety and portability checks

Checks run before a release, on a fresh clone so the working tree is not
rebuilt. Scripts:

| script | what |
|---|---|
| `run_matrix.sh` | gcc and clang at -O0/-O2/-O3, ASan+UBSan (both compilers), clang integer/bounds sanitizers, MSan, TSan; each runs `c/tests make check` against the frozen references; then valgrind memcheck on every test program |
| `run_static.sh` | cppcheck (warning, performance, portability; exhaustive) and gcc `-fanalyzer` over `c/src` and `c/wig` |
| `run_fuzz.sh`, `fuzz_rank1.c` | libFuzzer on `hexpr_rank1_pair` with the library instrumented (ASan+UBSan): status or in-range terms, at most `hexpr_rank1_max_terms`, the adjoint relation, the count-only call; then the line coverage of `csf_rank1.c` reached |
| `c/tests/test_oom.c` (in `make check`, Linux) | fails each allocation of ten public calls in turn (6706 failure points): out of memory reported, or the normal result; a normal run afterwards unchanged. Under the matrix's ASan/LSan and MSan builds this also checks the failure paths for leaks and invalid accesses |

## Results after the out-of-memory fix (release/v1.2.0, 2026-10-17)

The unchecked allocations found below were fixed (CHANGELOG, 1.2.0
"Fixed"). The matrix, valgrind and fuzzing ran on the fix (852d2f7); the
small follow-up that addressed the analyzer's remaining reports (inline
vector growth, a zeroed CSF work buffer, a dhash rehash that restores the
old table) was rechecked with `make check` under gcc and clang, ASan+UBSan
(both), MSan and `test_oom`. On the fixed code:

- **Matrix**: all 11 configurations pass `make check`, now including
  `test_oom`; no sanitizer report, no compiler warning; the frozen
  references are unchanged. valgrind: no failing program.
- **Fuzzing**: 25.8 million inputs in 300 s, no failure; coverage of
  `csf_rank1.c` as before.
- **Static analysis**: cppcheck reports nothing. gcc `-fanalyzer` reports
  one "leak" in `le_push` (`expr_tree.c`) on a path where `realloc`
  succeeds and its result is then taken to be NULL; the path is infeasible
  (the same code in `ct_push` / `ot_push` is not flagged) and is a false
  positive of gcc 14's realloc model.
- The bindings (Python, Julia, Ruby, Fortran) and hexpr2 dev/v1.1.0 (C gate
  and Python tests) pass against the fixed library.

## Results for 1.2.0 before the fix (release/v1.2.0 at c91cfb6, 2026-10-16)

Debian 13: gcc 14.2.0, clang 19.1.7, valgrind, cppcheck 2.17.1.

- **Matrix**: all 11 configurations pass `make check` (22 test programs and
  the Fortran-seam coverage script; outputs equal to the frozen references); no sanitizer report; no compiler warning
  with `-Wall -Wextra -pedantic`. The instrumentation was confirmed in the
  library and the test programs (`nm`).
- **valgrind**: no error and no definite or indirect leak in any of the 21
  test programs (including `test_pair_aggregate s_full_4o4e t_full_4o4e`).
- **Fuzzing**: 29.3 million inputs in 300 s, no failure. Coverage of
  `csf_rank1.c`: 97 % of lines; the only lines not reached are the error
  branch for `norb` above 256 (the harness keeps `norb` below 40).
- **Static analysis**: nothing in the code added in 1.2 (`csf_rank1.c`,
  `hexpr_lk_value_direct` in `expr_tree.c`, the new F77 wrappers). In the
  code of 1.0/1.1 both tools report unchecked allocations: a `malloc`,
  `calloc` or `realloc` whose NULL result is dereferenced, in `drt.c`,
  `expr.c`, `expr_tree.c` (`le_push`), `expr_internal.h`, `dhash.c` and
  `csf_pair_core.c`. They matter only when memory runs out; the calls and
  their output under normal conditions are not affected. gcc also flags the
  allocation size in `prune_unreachable` (`drt.c:657`, an `int` count whose
  range it cannot bound) and the `memcpy` from `ctx->csf_buf` in `csf_push`
  (`drt.c:1016`), whose buffer comes from an unchecked allocation.

Not done yet: Frama-C / CBMC on small core modules.
