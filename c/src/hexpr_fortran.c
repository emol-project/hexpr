/* hexpr_fortran.c -- F77 underscore-suffixed wrappers for libhexpr public API.
 * Not declared in hexpr.h; F77 callers link by name. Prototypes live in the
 * internal hexpr_fortran.h (included below) so the definitions are clean
 * under -Wmissing-prototypes. */
/**
 * @file hexpr_fortran.c
 * @brief F77 underscore-suffixed binding layer over the public @c hexpr_*() API.
 *
 * Each wrapper marshals Fortran-passed pointers (Fortran passes everything by
 * reference) into a single public @c hexpr_*() call and writes any results back
 * through caller-supplied out-pointers. There is no algorithm here -- only
 * dereference, forward, and write-back.
 *
 * @note Status convention: the @c void wrappers that can fail return the public
 *       ::hexpr_status_t through their @c *out_status argument, cast to @c int.
 *       This is stated once here and not repeated on every status wrapper.
 * @note Index convention: see the per-pair/block/subspace section note below --
 *       index arguments are forwarded WITHOUT base conversion.
 */
#include "../include/hexpr.h"
#include "hexpr_fortran.h"
#include <stdint.h>
#include <stdio.h>

/* ===== Fortran helpers ======================================= */

/* Returns the C stdout FILE* so the F90 *_show* wrappers can pass it to
 * hexpr_*_show*(..., FILE *).  Fortran cannot reliably alias the libc
 * "stdout" data symbol directly (copy relocation yields a zero-valued
 * local definition), so it calls this accessor instead. */
/** @brief Return the C @c stdout @c FILE* for the F90 @c *_show* wrappers to
 *  pass into @c hexpr_*_show*(..., FILE*). @return @c stdout as @c void*. */
void *hexpr_fortran_stdout(void) { return stdout; }

/* ===== Versioning ============================================ */

/** @brief Fortran binding for @ref hexpr_version_string. @return static version
 *  string (caller must NOT free). */
const char *hexpr_version_string_(void) { return hexpr_version_string(); }
/** @brief Fortran binding for @ref hexpr_version_major. @param[out] out major version. */
void hexpr_version_major_(int *out) { *out = hexpr_version_major(); }
/** @brief Fortran binding for @ref hexpr_version_minor. @param[out] out minor version. */
void hexpr_version_minor_(int *out) { *out = hexpr_version_minor(); }
/** @brief Fortran binding for @ref hexpr_version_patch. @param[out] out patch version. */
void hexpr_version_patch_(int *out) { *out = hexpr_version_patch(); }

/* ===== Status / error ======================================== */

/** @brief Fortran binding for @ref hexpr_last_error. @return the calling thread's
 *  last error string (thread-local, library-owned; caller must NOT free). */
const char *hexpr_last_error_(void) { return hexpr_last_error(); }

/* ===== Lifecycle ============================================= */

/** @brief Fortran binding for @ref hexpr_init. @param[out] out_status the
 *  ::hexpr_status_t result as @c int. */
void hexpr_init_(int *out_status)
{
    *out_status = (int)hexpr_init();
}

/** @brief Fortran binding for @ref hexpr_shutdown (no status). */
void hexpr_shutdown_(void) { hexpr_shutdown(); }

/* ===== DRT =================================================== */

/** @brief Fortran binding for @ref hexpr_drt_create. @param kind,nve,spin_x2,ncore,nval,next
 *  the by-reference scalar parameters. @return the new DRT handle (release with
 *  @ref hexpr_drt_destroy_), or NULL on failure. */
hexpr_drt_t *hexpr_drt_create_(int *kind, int *nve, int *spin_x2,
                                int *ncore, int *nval, int *next)
{
    return hexpr_drt_create((hexpr_drt_kind_t)*kind,
                            *nve, *spin_x2, *ncore, *nval, *next);
}

/** @brief Fortran binding for @ref hexpr_drt_from_csfset. @param norb,csfs,n_csfs
 *  the by-reference inputs. @return the new sub-DRT handle, or NULL on failure. */
hexpr_drt_t *hexpr_drt_from_csfset_(const int *norb,
                                     const unsigned char *csfs,
                                     const int *n_csfs)
{
    return hexpr_drt_from_csfset(*norb, csfs, *n_csfs);
}

/** @brief Fortran binding for @ref hexpr_drt_destroy. @param drt handle-by-reference;
 *  the pointed-to DRT is freed. */
void hexpr_drt_destroy_(hexpr_drt_t **drt)
{
    hexpr_drt_destroy(*drt);
}

/** @brief Fortran binding for @ref hexpr_drt_norb. @param drt DRT-by-reference.
 *  @param[out] out orbital count. */
void hexpr_drt_norb_(const hexpr_drt_t **drt, int *out)
{ *out = hexpr_drt_norb(*drt); }

/** @brief Fortran binding for @ref hexpr_drt_nnodes. @param drt DRT-by-reference.
 *  @param[out] out node count. */
void hexpr_drt_nnodes_(const hexpr_drt_t **drt, int *out)
{ *out = hexpr_drt_nnodes(*drt); }

/** @brief Fortran binding for @ref hexpr_drt_ncsfs. @param drt DRT-by-reference.
 *  @param[out] out user-facing CSF count. */
void hexpr_drt_ncsfs_(const hexpr_drt_t **drt, int *out)
{ *out = hexpr_drt_ncsfs(*drt); }

/** @brief Fortran binding for @ref hexpr_drt_walk_ncsfs. @param drt DRT-by-reference.
 *  @param[out] out canonical walk CSF count. */
void hexpr_drt_walk_ncsfs_(const hexpr_drt_t **drt, int *out)
{ *out = hexpr_drt_walk_ncsfs(*drt); }

/** @brief Fortran binding for @ref hexpr_drt_kind. @param drt DRT-by-reference.
 *  @param[out] out kind enum as @c int. */
void hexpr_drt_kind_(const hexpr_drt_t **drt, int *out)
{ *out = (int)hexpr_drt_kind(*drt); }

/** @brief Fortran binding for @ref hexpr_drt_ncore. @param drt DRT-by-reference.
 *  @param[out] out core orbital count. */
void hexpr_drt_ncore_(const hexpr_drt_t **drt, int *out)
{ *out = hexpr_drt_ncore(*drt); }

/** @brief Fortran binding for @ref hexpr_drt_nval. @param drt DRT-by-reference.
 *  @param[out] out valence orbital count. */
void hexpr_drt_nval_(const hexpr_drt_t **drt, int *out)
{ *out = hexpr_drt_nval(*drt); }

/** @brief Fortran binding for @ref hexpr_drt_next. @param drt DRT-by-reference.
 *  @param[out] out external orbital count. */
void hexpr_drt_next_(const hexpr_drt_t **drt, int *out)
{ *out = hexpr_drt_next(*drt); }

/** @brief Fortran binding for @ref hexpr_drt_nve_state. @param drt DRT-by-reference.
 *  @param[out] out target-state valence-electron count. */
void hexpr_drt_nve_state_(const hexpr_drt_t **drt, int *out)
{ *out = hexpr_drt_nve_state(*drt); }

/** @brief Fortran binding for @ref hexpr_drt_spin_x2. @param drt DRT-by-reference.
 *  @param[out] out target spin times two. */
void hexpr_drt_spin_x2_(const hexpr_drt_t **drt, int *out)
{ *out = hexpr_drt_spin_x2(*drt); }

/* CSF buffer: copies ncsfs*norb bytes into caller-supplied Fortran array */
/**
 * @brief Fortran binding for @ref hexpr_drt_csfs: copy the user-facing CSF buffer
 *        into a caller-supplied array.
 * @param drt DRT-by-reference.
 * @param[out] buf destination array.
 * @note Size contract: @p buf must hold at least @c ncsfs*norb bytes
 *       (@ref hexpr_drt_ncsfs_ x @ref hexpr_drt_norb_); exactly that many bytes
 *       are copied. A short buffer is a caller error.
 */
void hexpr_drt_csfs_(const hexpr_drt_t **drt, unsigned char *buf)
{
    const hexpr_drt_t *d = *drt;
    int n = hexpr_drt_ncsfs(d) * hexpr_drt_norb(d);
    const unsigned char *p = hexpr_drt_csfs(d);
    for (int i = 0; i < n; i++) buf[i] = p[i];
}

/* Canonical GUGA-walk CSF buffer: copies walk_ncsfs*norb bytes into the
 * caller-supplied Fortran array (CSFSET superset/canonical view). */
/**
 * @brief Fortran binding for @ref hexpr_drt_walk_csfs: copy the canonical
 *        GUGA-walk CSF buffer into a caller-supplied array.
 * @param drt DRT-by-reference.
 * @param[out] buf destination array.
 * @note Size contract: @p buf must hold at least @c walk_ncsfs*norb bytes
 *       (@ref hexpr_drt_walk_ncsfs_ x @ref hexpr_drt_norb_); exactly that many
 *       bytes are copied. A short buffer is a caller error.
 */
void hexpr_drt_walk_csfs_(const hexpr_drt_t **drt, unsigned char *buf)
{
    const hexpr_drt_t *d = *drt;
    int n = hexpr_drt_walk_ncsfs(d) * hexpr_drt_norb(d);
    const unsigned char *p = hexpr_drt_walk_csfs(d);
    for (int i = 0; i < n; i++) buf[i] = p[i];
}

/** @brief Fortran binding for @ref hexpr_drt_node_info. @param drt DRT-by-reference.
 *  @param node_index node index (by reference). @param[out] out_orb,out_nve,out_is
 *  the node code. @param[out] out_status status as @c int. */
void hexpr_drt_node_info_(const hexpr_drt_t **drt, int *node_index,
                           int *out_orb, int *out_nve, int *out_is,
                           int *out_status)
{
    *out_status = (int)hexpr_drt_node_info(*drt, *node_index,
                                            out_orb, out_nve, out_is);
}

/** @brief Fortran binding for @ref hexpr_drt_show, printing to C @c stdout.
 *  @param drt DRT-by-reference. @param[out] out_status status as @c int. */
void hexpr_drt_show_(const hexpr_drt_t **drt, int *out_status)
{
    *out_status = (int)hexpr_drt_show(*drt, stdout);
}

/* ===== Expression ============================================ */

/** @brief Fortran binding for @ref hexpr_expr_build. @param drt DRT-by-reference.
 *  @return the expression handle (release with @ref hexpr_expr_destroy_). */
hexpr_expr_t *hexpr_expr_build_(const hexpr_drt_t **drt)
{
    return hexpr_expr_build(*drt);
}

/** @brief Fortran binding for @ref hexpr_expr_destroy. @param expr handle-by-reference;
 *  the pointed-to expression is freed. */
void hexpr_expr_destroy_(hexpr_expr_t **expr)
{
    hexpr_expr_destroy(*expr);
}

/** @brief Fortran binding for @ref hexpr_expr_nterms_full. @param expr handle-by-reference.
 *  @param[out] out total term count. */
void hexpr_expr_nterms_full_(const hexpr_expr_t **expr, int64_t *out)
{ *out = hexpr_expr_nterms_full(*expr); }

/** @brief Fortran binding for @ref hexpr_expr_nterms_one. @param expr handle-by-reference.
 *  @param[out] out one-electron term count. */
void hexpr_expr_nterms_one_(const hexpr_expr_t **expr, int64_t *out)
{ *out = hexpr_expr_nterms_one(*expr); }

/** @brief Fortran binding for @ref hexpr_expr_read. @param expr handle-by-reference.
 *  @param which 0=full / 1=one-electron (by reference). @param[out] out_ij,out_pqrs,out_coef
 *  parallel term columns. @param[out] out_status status as @c int. */
void hexpr_expr_read_(const hexpr_expr_t **expr, int *which,
                      int32_t *out_ij, int32_t *out_pqrs, double *out_coef,
                      int *out_status)
{
    *out_status = (int)hexpr_expr_read(*expr, (hexpr_expr_which_t)*which,
                                        out_ij, out_pqrs, out_coef);
}

/** @brief Fortran binding for @ref hexpr_expr_show_full, printing to C @c stdout.
 *  @param expr handle-by-reference. @param[out] out_status status as @c int. */
void hexpr_expr_show_full_(const hexpr_expr_t **expr, int *out_status)
{ *out_status = (int)hexpr_expr_show_full(*expr, stdout); }

/** @brief Fortran binding for @ref hexpr_expr_show_one, printing to C @c stdout.
 *  @param expr handle-by-reference. @param[out] out_status status as @c int. */
void hexpr_expr_show_one_(const hexpr_expr_t **expr, int *out_status)
{ *out_status = (int)hexpr_expr_show_one(*expr, stdout); }

/* ===== CSF lookup ============================================ */

/* step is an array of Fortran integer(c_int8_t) (signed char) step codes;
 * cast to unsigned char for the public API. */
/**
 * @brief Fortran binding for @ref hexpr_csf_index.
 * @param drt DRT-by-reference. @param step step-code array. @param[out] out_index
 *        the resulting index. @param[out] out_status status as @c int.
 * @note The Fortran @c integer(c_int8_t) step array is @c signed @c char here and
 *       is reinterpret-cast to the public API's @c unsigned @c char*.
 */
void hexpr_csf_index_(const hexpr_drt_t **drt, const signed char *step,
                       int *out_index, int *out_status)
{
    *out_status = (int)hexpr_csf_index(*drt, (const unsigned char *)step,
                                        out_index);
}

/* Writes norb step-code bytes into the caller-supplied out_step buffer. */
/**
 * @brief Fortran binding for @ref hexpr_csf_steps.
 * @param drt DRT-by-reference. @param index CSF index (by reference).
 * @param[out] out_step receives @c norb step-code bytes. @param[out] out_status
 *        status as @c int.
 * @note The @c out_step array is Fortran @c integer(c_int8_t) (@c signed @c char
 *       here) and is reinterpret-cast to the public API's @c unsigned @c char*.
 * @note @c *index is forwarded without base conversion (see the
 *       per-pair/block/subspace section note below).
 */
void hexpr_csf_steps_(const hexpr_drt_t **drt, const int *index,
                       signed char *out_step, int *out_status)
{
    *out_status = (int)hexpr_csf_steps(*drt, *index,
                                        (unsigned char *)out_step);
}

/* Convention: Fortran passes the CHARACTER buffer (str) and an explicit
 * orbital count (norb); the hidden trailing string-length argument that
 * Fortran appends after all explicit arguments is ignored.  The first
 * *norb characters of str are copied into a NUL-terminated local buffer
 * before calling hexpr_step_parse. */
/**
 * @brief Fortran binding for @ref hexpr_step_parse: parse an "eudf" CHARACTER
 *        buffer into step codes.
 * @param str the Fortran CHARACTER buffer. @param norb orbital count (by reference).
 * @param[out] out_step receives the parsed step codes. @param[out] out_status
 *        status as @c int.
 * @note The only wrapper with logic of its own: it clamps @c *norb to @c [0,255],
 *       copies that many characters from @p str into a fixed 256-byte local
 *       buffer, and NUL-terminates it before calling @ref hexpr_step_parse.
 * @note The hidden trailing Fortran string-length argument (appended after all
 *       explicit arguments) is intentionally ignored; the caller passes @p norb
 *       explicitly. (Formalises the convention block comment above.)
 */
void hexpr_step_parse_(const char *str, const int *norb,
                        signed char *out_step, int *out_status)
{
    char buf[256];
    int n = *norb;
    if (n < 0)   n = 0;
    if (n > 255) n = 255;
    for (int i = 0; i < n; i++) buf[i] = str[i];
    buf[n] = '\0';
    *out_status = (int)hexpr_step_parse(buf, *norb,
                                         (unsigned char *)out_step);
}

/* ===== CSF -> determinant expansion ========================== */

/**
 * @brief Fortran binding for @ref hexpr_csf_ndets.
 * @param drt DRT-by-reference. @param index CSF index (by reference).
 * @param ms_x2 2*Ms (by reference). @param[out] out_ndets determinant count.
 * @param[out] out_status status as @c int.
 * @note @c *index is forwarded without base conversion (see the
 *       per-pair/block/subspace section note below).
 */
void hexpr_csf_ndets_(const hexpr_drt_t **drt, const int *index,
                       const int *ms_x2, int *out_ndets, int *out_status)
{
    *out_status = (int)hexpr_csf_ndets(*drt, *index, *ms_x2, out_ndets);
}

/**
 * @brief Fortran binding for @ref hexpr_csf_to_dets.
 * @param drt DRT-by-reference. @param index CSF index (by reference).
 * @param ms_x2 2*Ms (by reference). @param[out] out_count determinants written.
 * @param[out] out_alpha,out_beta occupation bit masks.
 * @param[out] out_coef expansion coefficients. @param[out] out_status status
 *        as @c int.
 * @note Size contract: @p out_alpha and @p out_beta must each hold at least
 *       @c ndets*nwords elements and @p out_coef at least @c ndets, where
 *       @c ndets comes from @ref hexpr_csf_ndets_ and
 *       @c nwords = (norb+63)/64. A short buffer is a caller error.
 * @note Fortran has no unsigned integer type: the masks arrive as
 *       @c integer*8 / @c integer(c_int64_t), so a determinant occupying
 *       spatial orbital 63 sets the sign bit and the value reads as
 *       negative. Test occupation with @c BTEST, never with a magnitude
 *       comparison.
 * @note @c *index is forwarded without base conversion (see the
 *       per-pair/block/subspace section note below).
 */
void hexpr_csf_to_dets_(const hexpr_drt_t **drt, const int *index,
                         const int *ms_x2, int *out_count,
                         uint64_t *out_alpha, uint64_t *out_beta,
                         double *out_coef, int *out_status)
{
    *out_status = (int)hexpr_csf_to_dets(*drt, *index, *ms_x2, out_count,
                                          out_alpha, out_beta, out_coef);
}

/* ===== Per-pair / block / subspace =========================== */

/**
 * @name Per-pair / block / subspace bindings
 * @note Base-index contract: this binding performs NO base conversion on any
 *       index argument. The @c bra / @c ket / @c index / @c indices values are
 *       forwarded to the public API unchanged. (Statement of absence: no
 *       @c +1/-1 index adjustment occurs in this file; the @c +1 inside the
 *       buffer-copy loops above is a length bound, not an index.) This scopes
 *       the index-taking wrappers below (count/eval families and
 *       @ref hexpr_csf_steps_); it makes no claim about any F90 layer.
 * @{
 */

/** @brief Fortran binding for @ref hexpr_pair_count.
 *  @param drt DRT-by-reference. @param bra,ket,which by-reference scalars.
 *  @param[out] out_count term count. @param[out] out_status status as @c int.
 *  @note index forwarded unchanged, no base conversion (see section note). */
void hexpr_pair_count_(const hexpr_drt_t **drt,
                        const int *bra, const int *ket,
                        const int *which, int *out_count, int *out_status)
{
    *out_status = (int)hexpr_pair_count(*drt, *bra, *ket,
                                         (hexpr_expr_which_t)*which, out_count);
}

/** @brief Fortran binding for @ref hexpr_pair_max_terms.
 *  @param drt DRT-by-reference. @param which by-reference scalar.
 *  @param[out] out_max term upper bound. @param[out] out_status status as @c int. */
void hexpr_pair_max_terms_(const hexpr_drt_t **drt,
                             const int *which, int *out_max, int *out_status)
{
    *out_status = (int)hexpr_pair_max_terms(*drt,
                                             (hexpr_expr_which_t)*which, out_max);
}

/** @brief Fortran binding for @ref hexpr_pair_eval.
 *  @param drt DRT-by-reference. @param bra,ket,which by-reference scalars.
 *  @param[out] out_count terms written. @param[out] out_ij,out_pqrs,out_coef
 *  parallel columns. @param[out] out_status status as @c int.
 *  @note index forwarded unchanged, no base conversion (see section note). */
void hexpr_pair_eval_(const hexpr_drt_t **drt,
                       const int *bra, const int *ket, const int *which,
                       int *out_count,
                       int32_t *out_ij, int32_t *out_pqrs, double *out_coef,
                       int *out_status)
{
    *out_status = (int)hexpr_pair_eval(*drt, *bra, *ket,
                                        (hexpr_expr_which_t)*which,
                                        out_count, out_ij, out_pqrs, out_coef);
}

/** @brief Fortran binding for @ref hexpr_block_count.
 *  @param drt DRT-by-reference. @param bra_indices,n_bra,ket_indices,n_ket,which
 *  by-reference inputs. @param[out] out_count term count. @param[out] out_status
 *  status as @c int.
 *  @note index forwarded unchanged, no base conversion (see section note). */
void hexpr_block_count_(const hexpr_drt_t **drt,
                         const int *bra_indices, const int *n_bra,
                         const int *ket_indices, const int *n_ket,
                         const int *which, int *out_count, int *out_status)
{
    *out_status = (int)hexpr_block_count(*drt, bra_indices, *n_bra,
                                          ket_indices, *n_ket,
                                          (hexpr_expr_which_t)*which, out_count);
}

/** @brief Fortran binding for @ref hexpr_block_max_terms.
 *  @param drt DRT-by-reference. @param n_bra,n_ket,which by-reference scalars.
 *  @param[out] out_max term upper bound. @param[out] out_status status as @c int. */
void hexpr_block_max_terms_(const hexpr_drt_t **drt,
                              const int *n_bra, const int *n_ket,
                              const int *which, int *out_max, int *out_status)
{
    *out_status = (int)hexpr_block_max_terms(*drt, *n_bra, *n_ket,
                                              (hexpr_expr_which_t)*which, out_max);
}

/** @brief Fortran binding for @ref hexpr_block_eval.
 *  @param drt DRT-by-reference. @param bra_indices,n_bra,ket_indices,n_ket,which
 *  by-reference inputs. @param[out] out_count terms written. @param[out]
 *  out_bra_pos,out_ket_pos set positions. @param[out] out_ij,out_pqrs,out_coef
 *  parallel columns. @param[out] out_status status as @c int.
 *  @note index forwarded unchanged, no base conversion (see section note). */
void hexpr_block_eval_(const hexpr_drt_t **drt,
                        const int *bra_indices, const int *n_bra,
                        const int *ket_indices, const int *n_ket,
                        const int *which, int *out_count,
                        int32_t *out_bra_pos, int32_t *out_ket_pos,
                        int32_t *out_ij, int32_t *out_pqrs, double *out_coef,
                        int *out_status)
{
    *out_status = (int)hexpr_block_eval(*drt, bra_indices, *n_bra,
                                         ket_indices, *n_ket,
                                         (hexpr_expr_which_t)*which, out_count,
                                         out_bra_pos, out_ket_pos,
                                         out_ij, out_pqrs, out_coef);
}

/** @brief Fortran binding for @ref hexpr_subspace_count.
 *  @param drt DRT-by-reference. @param indices,n,which by-reference inputs.
 *  @param[out] out_count term count. @param[out] out_status status as @c int.
 *  @note index forwarded unchanged, no base conversion (see section note). */
void hexpr_subspace_count_(const hexpr_drt_t **drt,
                             const int *indices, const int *n,
                             const int *which, int *out_count, int *out_status)
{
    *out_status = (int)hexpr_subspace_count(*drt, indices, *n,
                                             (hexpr_expr_which_t)*which, out_count);
}

/** @brief Fortran binding for @ref hexpr_subspace_max_terms.
 *  @param drt DRT-by-reference. @param n,which by-reference scalars.
 *  @param[out] out_max term upper bound. @param[out] out_status status as @c int. */
void hexpr_subspace_max_terms_(const hexpr_drt_t **drt,
                                 const int *n, const int *which,
                                 int *out_max, int *out_status)
{
    *out_status = (int)hexpr_subspace_max_terms(*drt, *n,
                                                 (hexpr_expr_which_t)*which, out_max);
}

/** @brief Fortran binding for @ref hexpr_subspace_eval.
 *  @param drt DRT-by-reference. @param indices,n,which by-reference inputs.
 *  @param[out] out_count terms written. @param[out] out_bra_pos,out_ket_pos set
 *  positions. @param[out] out_ij,out_pqrs,out_coef parallel columns. @param[out]
 *  out_status status as @c int.
 *  @note index forwarded unchanged, no base conversion (see section note). */
void hexpr_subspace_eval_(const hexpr_drt_t **drt,
                           const int *indices, const int *n,
                           const int *which, int *out_count,
                           int32_t *out_bra_pos, int32_t *out_ket_pos,
                           int32_t *out_ij, int32_t *out_pqrs, double *out_coef,
                           int *out_status)
{
    *out_status = (int)hexpr_subspace_eval(*drt, indices, *n,
                                            (hexpr_expr_which_t)*which, out_count,
                                            out_bra_pos, out_ket_pos,
                                            out_ij, out_pqrs, out_coef);
}
/** @} */
