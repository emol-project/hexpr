/* csfset_oracle_util.h -- shared constants for the CSFSET pair-eval oracle.
 *
 * Used by gen_csfset_oracle.c (generator) and test_csfset_pair.c (consumer).
 *
 * This header used to carry csfset_canon_index(), which mapped a CSFSET input
 * row to its canonical (walk) index by scanning hexpr_drt_walk_csfs(), because
 * the eval API addressed the canonical superset while bra_pos/ket_pos reported
 * positions in the array we passed. Both callers had to apply it, and the same
 * mapping was written three more times in the Python, Julia and Ruby bindings.
 * The eval API now takes INPUT-ORDER indices and does the mapping itself, so
 * the index array is simply 0..n-1 and the helper is gone.
 */
#ifndef CSFSET_ORACLE_UTIL_H_
#define CSFSET_ORACLE_UTIL_H_

#include "hexpr.h"

/* Directory of committed oracle files, relative to c/tests/ run dir. */
#define CSFSET_ORACLE_DIR "../../reference_csfset"

/* The six frozen oracle cases (see docs/SPEC_CSFSET_ORACLE.md). */
static const char *const CSFSET_ORACLE_CASES[] = {
    "O1", "O2", "O5", "OF", "OFs", "OS"
};
enum { CSFSET_ORACLE_NCASES =
    (int)(sizeof(CSFSET_ORACLE_CASES) / sizeof(CSFSET_ORACLE_CASES[0])) };

#endif /* CSFSET_ORACLE_UTIL_H_ */
