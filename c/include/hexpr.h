#ifndef HEXPR_H_
#define HEXPR_H_

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Everything declared below is the installed ABI. The library is built with
 * -fvisibility=hidden, so this region is what makes these symbols -- and
 * only these -- visible in libhexpr.so. System headers must stay outside
 * it: a visibility push applies to nested includes as well. */
#pragma GCC visibility push(default)

/* ===== Versioning ============================================ */
#define HEXPR_VERSION_MAJOR 1
#define HEXPR_VERSION_MINOR 0
#define HEXPR_VERSION_PATCH 0

const char *hexpr_version_string(void);
int         hexpr_version_major (void);
int         hexpr_version_minor (void);
int         hexpr_version_patch (void);

/* ===== Status codes ========================================== */
typedef enum {
    HEXPR_OK              = 0,
    HEXPR_ERR_INVALID_ARG = 1,
    HEXPR_ERR_OOM         = 2,
    HEXPR_ERR_IO          = 3,
    HEXPR_ERR_RANGE       = 4,
    HEXPR_ERR_NOT_BUILT   = 5,
    HEXPR_ERR_INTERNAL    = 99
} hexpr_status_t;

/* Thread-local. Pointer valid until next hexpr call from same thread. */
const char *hexpr_last_error(void);

/* ===== Library lifecycle ===================================== */

/* Idempotent. Calling again returns HEXPR_OK without doing anything. */
hexpr_status_t hexpr_init(void);

void hexpr_shutdown(void);

/* ===== DRT =================================================== */

typedef struct hexpr_drt hexpr_drt_t;

typedef enum {
    HEXPR_DRT_SOCI   = 0,
    HEXPR_DRT_FOCI   = 1,
    HEXPR_DRT_CSFSET = 2
} hexpr_drt_kind_t;

hexpr_drt_t *hexpr_drt_create(hexpr_drt_kind_t kind,
                              int nve,
                              int spin_x2,
                              int ncore,
                              int nval,
                              int next);

void hexpr_drt_destroy(hexpr_drt_t *drt);

/* Build a sub-DRT from an explicit CSF set. Each CSF must be length norb and
 * use step bytes 0..3 (HEXPR_STEP_E/U/D/F). All CSFs must share a common top
 * node. csfs is a packed buffer of n_csfs * norb bytes.
 * Returns NULL on validation failure (hexpr_last_error() is set).
 *
 * User-facing CSF view (input order): for the returned CSFSET sub-DRT the
 * user-facing accessors speak the INPUT set in INPUT order -- hexpr_drt_csfs()
 * returns exactly these n_csfs rows in the supplied order, hexpr_drt_ncsfs()
 * == n_csfs, hexpr_csf_steps(i) is input row i, and hexpr_csf_index() returns
 * the input position. Internally the sub-DRT re-derives a canonical GUGA-walk
 * CSF list (e<u<d<f) that may be a SUPERSET (more CSFs than supplied, when
 * input CSFs share nodes); that canonical view is available via
 * hexpr_drt_walk_csfs() / hexpr_drt_walk_ncsfs(). The bra_pos/ket_pos indices
 * from hexpr_block_eval / hexpr_subspace_eval index the input array and so
 * agree with hexpr_drt_csfs() row order. */
hexpr_drt_t *hexpr_drt_from_csfset(int norb,
                                    const unsigned char *csfs,
                                    int n_csfs);

/* Orbital/CSF counts */
int  hexpr_drt_norb     (const hexpr_drt_t *);
int  hexpr_drt_nnodes   (const hexpr_drt_t *);
int  hexpr_drt_ncsfs    (const hexpr_drt_t *);
hexpr_drt_kind_t hexpr_drt_kind(const hexpr_drt_t *);

/* MO classification (needed by FFI users to reconstruct full DRT metadata) */
int  hexpr_drt_ncore    (const hexpr_drt_t *);
int  hexpr_drt_nval     (const hexpr_drt_t *);
int  hexpr_drt_next     (const hexpr_drt_t *);

/* State parameters */
int  hexpr_drt_nve_state (const hexpr_drt_t *);  /* number of valence electrons */
int  hexpr_drt_spin_x2   (const hexpr_drt_t *);  /* 2*S (spin quantum number * 2) */

/* CSF list as packed buffer of step codes (E/U/D/F = 0/1/2/3). */
typedef enum {
    HEXPR_STEP_E = 0,
    HEXPR_STEP_U = 1,
    HEXPR_STEP_D = 2,
    HEXPR_STEP_F = 3
} hexpr_step_t;

/* User-facing CSF list: pointer to ncsfs * norb bytes, owned by DRT. DO NOT
 * FREE. For a CSFSET sub-DRT these are the n_csfs input rows in input order
 * (ncsfs == hexpr_drt_ncsfs()); for SOCI/FOCI they are the parametric CSFs.
 * The canonical (possibly superset) view is hexpr_drt_walk_csfs(). */
const unsigned char *hexpr_drt_csfs(const hexpr_drt_t *);

/* Canonical GUGA-walk CSF list and count, for all kinds. For a CSFSET
 * sub-DRT this is the internal view, which may be a SUPERSET of the input
 * (hexpr_drt_walk_ncsfs() >= hexpr_drt_ncsfs()); for SOCI/FOCI it equals
 * hexpr_drt_csfs() / hexpr_drt_ncsfs(). Buffer owned by DRT; DO NOT FREE. */
const unsigned char *hexpr_drt_walk_csfs(const hexpr_drt_t *);
int                  hexpr_drt_walk_ncsfs(const hexpr_drt_t *);

/* ===== CSF lookup ============================================ */

/* Step codes -> CSF index. step is length-norb unsigned char array; each
 * byte is 0..3 (matching hexpr_step_t values). For a CSFSET sub-DRT the
 * returned index is the INPUT position (0..ncsfs-1, input order); a valid
 * step path that is NOT one of the input CSFs (e.g. a canonical superset
 * member) returns HEXPR_ERR_RANGE. For SOCI/FOCI it is the canonical DRT
 * row index. */
hexpr_status_t hexpr_csf_index(const hexpr_drt_t *drt,
                                const unsigned char *step,
                                int *out_index);

/* CSF index -> step codes. out_step must point to a length-norb unsigned
 * char buffer. For a CSFSET sub-DRT, index is the INPUT position
 * (0..ncsfs-1) and the result is that input row; for SOCI/FOCI it is the
 * canonical row. */
hexpr_status_t hexpr_csf_steps(const hexpr_drt_t *drt,
                                int index,
                                unsigned char *out_step);

/* Parse step string ("eudfee...") into byte array. Case-insensitive;
 * non-{e,u,d,f} characters cause HEXPR_ERR_INVALID_ARG.
 * out_step must have length >= norb. */
hexpr_status_t hexpr_step_parse(const char *str, int norb,
                                 unsigned char *out_step);

/* ===== CSF -> determinant expansion ========================== */

/* Expand a CSF into its Slater-determinant components:
 *
 *     |CSF_I> = sum_k c_k |det_k>
 *
 * Conventions, all fixed here and none inferred from a passing test.
 * src/csf/csf2det.c states them at greater length and derives them.
 *
 *   Spin coupling: the genealogical coupling defined by the DRT walk.
 *     A CSF expressed in another coupling scheme will not reproduce
 *     these coefficients.
 *
 *   Clebsch-Gordan phase: Condon-Shortley.
 *
 *   Operator order: ALL ALPHA FIRST in ascending spatial orbital, THEN
 *     ALL BETA in ascending spatial orbital.
 *
 *       |det> = a^dag_{p1,alpha} ... a^dag_{pm,alpha}
 *               a^dag_{q1,beta}  ... a^dag_{qn,beta} |vac>,
 *       p1 < ... < pm,  q1 < ... < qn.
 *
 *     This is the block convention used by most determinant CI codes
 *     and by most second-quantized qubit mappings. An
 *     interleaved-by-orbital convention differs from it by a
 *     determinant-dependent permutation sign; csf2det.c gives the rule.
 *
 *   Ms: passed as 2*Ms in ms_x2. Any |Ms| <= S with the same parity as
 *     2*S is accepted, negative values included. Ms is NOT assumed to
 *     equal S: a comparison against a determinant program run at
 *     nelec=(na,nb) fixes Ms = (na-nb)/2 on that side, and that side's
 *     choice is not ours to make.
 *
 *   Normalization: sum_k c_k^2 == 1 for every CSF, at every Ms.
 *
 *   Zero coefficients are emitted, not dropped, so the count is exact
 *   and no coefficient threshold enters the API.
 *
 *   Determinant order: lexicographically ascending in the set of
 *     open-shell positions (walk order, 0-based) that carry alpha.
 */

/* Number of determinants in the expansion of CSF `index` at `ms_x2`.
 * Exact, not an upper bound: it is C(n_open, n_alpha_open) and follows
 * from the step vector and ms_x2 alone.
 *
 * `index` follows the same convention as hexpr_csf_steps(): the INPUT
 * position for a CSFSET sub-DRT, the canonical row for SOCI/FOCI.
 *
 * Returns a status rather than the count directly, unlike the plain
 * accessors above, because callers size allocations with this value and
 * a silently wrong count is an out-of-bounds read at the call site. */
hexpr_status_t hexpr_csf_ndets(const hexpr_drt_t *drt,
                               int index,
                               int ms_x2,
                               int *out_ndets);

/* Expand CSF `index` into determinants at `ms_x2`.
 *
 * Let ndets = hexpr_csf_ndets(drt, index, ms_x2) and
 *     nwords = (hexpr_drt_norb(drt) + 63) / 64.
 *
 * out_alpha and out_beta must each point to a buffer of ndets * nwords
 * uint64_t; out_coef to a buffer of ndets double. Determinant k occupies
 * out_alpha[k*nwords .. k*nwords + nwords - 1], and likewise for beta.
 * Within a determinant, bit p of word w is spatial orbital 64*w + p.
 *
 * out_count receives the number of determinants written, which equals
 * hexpr_csf_ndets() for the same arguments. */
hexpr_status_t hexpr_csf_to_dets(const hexpr_drt_t *drt,
                                 int index,
                                 int ms_x2,
                                 int *out_count,
                                 uint64_t *out_alpha,
                                 uint64_t *out_beta,
                                 double   *out_coef);

hexpr_status_t hexpr_drt_node_info(const hexpr_drt_t *, int node_index,
                                   int *out_orb,
                                   int *out_nve,
                                   int *out_is);

/* Ruby text format for regression. Deprecated; for legacy regression only. */
hexpr_status_t hexpr_drt_show(const hexpr_drt_t *, FILE *out);

/* ===== Energy expression ===================================== */

/* --- Output index convention (ALL binary read APIs below: hexpr_expr_read,
 *     hexpr_pair_eval, hexpr_block_eval, hexpr_subspace_eval) ---
 *
 * out_ij is a 0-based PACKED triangle index of the (bra,ket) CSF pair -- NOT a
 * row index. Decode it to 0-based (I,J), I >= J:
 *     I = (int)floor((sqrt(8.0 * ij + 1.0) - 1.0) / 2.0);
 *     J = ij - I * (I + 1) / 2;
 * then apply any language +1 at the array-indexing layer (e.g. Fortran/Julia
 * 1-based arrays). This 0-based binary convention matches the PySCF/Psi4/Fermi.jl
 * integration seams, which decode 0-based and add the +1 in-language.
 *
 * out_pqrs is a 0-based packed integral index (seq_num): one-electron terms use
 * tri(p,q)-1 in [0, n_oel); two-electron terms use tri(tri(p,q),tri(r,s))-1+n_oel,
 * with n_oel = norb*(norb+1)/2 and tri(a,b) = max*(max-1)/2 + min (1-based p,q,r,s).
 *
 * The text dumps (hexpr_expr_show_full/_one) print ij 1-based (ij+1) to match the
 * Ruby reference corpus reference/<case>.expr.txt, but pqrs 0-based -- a deliberate
 * dump-only asymmetry. The binary APIs here emit BOTH ij and pqrs 0-based. */

typedef struct hexpr_expr hexpr_expr_t;

hexpr_expr_t *hexpr_expr_build(const hexpr_drt_t *);
void          hexpr_expr_destroy(hexpr_expr_t *);

int64_t hexpr_expr_nterms_full(const hexpr_expr_t *);
int64_t hexpr_expr_nterms_one (const hexpr_expr_t *);

typedef enum {
    HEXPR_EXPR_FULL = 0,
    HEXPR_EXPR_ONE  = 1
} hexpr_expr_which_t;

/* Bulk read. Caller allocates three parallel arrays of length
 * `nterms` (= hexpr_expr_nterms_full or _nterms_one).
 * Returns HEXPR_OK on success. */
hexpr_status_t hexpr_expr_read(const hexpr_expr_t *,
                               hexpr_expr_which_t which,
                               int32_t *out_ij,
                               int32_t *out_pqrs,
                               double  *out_coef);

/* Ruby text format for regression. Deprecated; for legacy regression only. */
hexpr_status_t hexpr_expr_show_full(const hexpr_expr_t *, FILE *out);
hexpr_status_t hexpr_expr_show_one (const hexpr_expr_t *, FILE *out);

/* ===== Per-pair expression evaluation ======================== */

/* The per-pair / block / subspace evaluators return one direction-selected
 * decomposition of H[I,J]; see docs/PAIR_API_SEMANTICS.md for the directional
 * semantics and how to obtain the complete decomposition. */

/* --- Index convention for the six calls below ---
 *
 * Every CSF index passed to hexpr_pair_count/_eval, hexpr_block_count/_eval and
 * hexpr_subspace_count/_eval is an INPUT-ORDER index in [0, hexpr_drt_ncsfs()):
 * the same space as hexpr_drt_csfs(), hexpr_csf_index() and hexpr_csf_steps().
 * There is one index space in this API. For SOCI and FOCI it is also the
 * canonical walk order; for a CSFSET sub-DRT the canonical walk may reorder,
 * extend or dedup the input, and the translation is done inside the library.
 *
 * out_ij is the ONE exception, and it is not an identifier a caller needs here.
 * Per docs/PQRS_INDEX_SPEC.md (normative for v1.0.0) ij is triangle() of the two
 * Shavitt PATH WEIGHTS, so it stays in the canonical space; an input position is
 * not a path weight. To attribute a term to a pair, use out_bra_pos/out_ket_pos
 * from the block and subspace forms -- they are positions within the arrays you
 * passed. Do NOT decode out_ij into CSF indices and feed them back to these
 * calls or to hexpr_csf_steps(); for a reordered CSFSET sub-DRT they name
 * different CSFs. For SOCI/FOCI, and for hexpr_expr_read() over a full DRT,
 * decoding ij is exactly the row/column pair as documented above. */

/* Count expression terms for a single (bra, ket) CSF pair.
 * bra/ket are INPUT-ORDER CSF indices; see the note above. */
hexpr_status_t hexpr_pair_count(const hexpr_drt_t *drt, int bra, int ket,
                                 hexpr_expr_which_t which, int *out_count);

/* Evaluate a single (bra, ket) CSF pair.
 * Caller allocates out_ij/out_pqrs/out_coef with hexpr_pair_count elements
 * (any of these may be NULL if not needed).
 * *out_count is set to the number of terms written.
 * bra/ket are INPUT-ORDER CSF indices, in [0, hexpr_drt_ncsfs()); see the note
 * above. The emitted out_ij stays a canonical pair address. */
hexpr_status_t hexpr_pair_eval(const hexpr_drt_t *drt, int bra, int ket,
                                hexpr_expr_which_t which,
                                int *out_count,
                                int32_t *out_ij, int32_t *out_pqrs,
                                double *out_coef);

/* ===== Block evaluation (rectangular bra x ket) ============= */

/* Count expression terms over a rectangular bra x ket block.
 * bra_indices/ket_indices are INPUT-ORDER CSF indices; see the note above. */
hexpr_status_t hexpr_block_count(const hexpr_drt_t *drt,
                                  const int *bra_indices, int n_bra,
                                  const int *ket_indices, int n_ket,
                                  hexpr_expr_which_t which,
                                  int *out_count);

/* Evaluate a rectangular bra x ket block of CSF pairs.
 * bra_indices/ket_indices are INPUT-ORDER CSF indices, in
 * [0, hexpr_drt_ncsfs()); see the note above. The returned
 * out_bra_pos/out_ket_pos index the caller's bra_indices/ket_indices arrays, in
 * the order supplied, and are how a term is attributed to a pair. */
hexpr_status_t hexpr_block_eval(const hexpr_drt_t *drt,
                                 const int *bra_indices, int n_bra,
                                 const int *ket_indices, int n_ket,
                                 hexpr_expr_which_t which,
                                 int *out_count,
                                 int32_t *out_bra_pos,
                                 int32_t *out_ket_pos,
                                 int32_t *out_ij,
                                 int32_t *out_pqrs,
                                 double  *out_coef);

/* ===== Subspace evaluation (bra == ket) ===================== */

/* Count expression terms over a square subspace (bra == ket index set).
 * indices are INPUT-ORDER CSF indices; see the note above. */
hexpr_status_t hexpr_subspace_count(const hexpr_drt_t *drt,
                                     const int *indices, int n,
                                     hexpr_expr_which_t which,
                                     int *out_count);

/* Evaluate a square subspace (bra == ket index set) of CSF pairs.
 * indices are INPUT-ORDER CSF indices, in [0, hexpr_drt_ncsfs()); see the note
 * above. The returned out_bra_pos/out_ket_pos index the caller's indices array,
 * in the order supplied, and are how a term is attributed to a pair. */
hexpr_status_t hexpr_subspace_eval(const hexpr_drt_t *drt,
                                    const int *indices, int n,
                                    hexpr_expr_which_t which,
                                    int *out_count,
                                    int32_t *out_bra_pos,
                                    int32_t *out_ket_pos,
                                    int32_t *out_ij,
                                    int32_t *out_pqrs,
                                    double  *out_coef);

/* ===== Max-terms upper bounds ================================= */

/* Returns an upper bound on the number of output terms a single
 * hexpr_pair_eval call can produce. Useful for sizing a static
 * buffer (Fortran 77 cannot use ALLOCATE). The actual number of
 * terms is typically much smaller. */
hexpr_status_t hexpr_pair_max_terms(const hexpr_drt_t *drt,
                                     hexpr_expr_which_t which,
                                     int *out_max);

/* Upper bound for hexpr_block_eval. Returns HEXPR_ERR_RANGE if
 * the upper bound exceeds INT_MAX (32-bit signed overflow). */
hexpr_status_t hexpr_block_max_terms(const hexpr_drt_t *drt,
                                      int n_bra, int n_ket,
                                      hexpr_expr_which_t which,
                                      int *out_max);

/* Upper bound for hexpr_subspace_eval. Returns HEXPR_ERR_RANGE if
 * the upper bound exceeds INT_MAX. */
hexpr_status_t hexpr_subspace_max_terms(const hexpr_drt_t *drt,
                                         int n,
                                         hexpr_expr_which_t which,
                                         int *out_max);

#pragma GCC visibility pop

#ifdef __cplusplus
}
#endif

#endif /* HEXPR_H_ */
