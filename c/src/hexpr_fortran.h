/* hexpr_fortran.h -- prototypes for the F77 underscore-suffixed wrappers
 * defined in hexpr_fortran.c.
 *
 * These symbols are NOT part of the public C API (hexpr.h); Fortran callers
 * link against them by name via their own iso_c_binding interfaces. This
 * header exists only so the wrapper definitions have visible prototypes
 * (prototype hygiene -- builds cleanly under -Wmissing-prototypes) and so the
 * ABI is documented in one place. It is internal: do not install it. */
#ifndef HEXPR_FORTRAN_H_
#define HEXPR_FORTRAN_H_

#include "hexpr.h"
#include <stdint.h>

/* The wrappers are not public C API, but Fortran callers bind to them by
 * name in libhexpr.so, so they must stay visible under -fvisibility=hidden.
 * Same rule as hexpr.h: no system header inside the region. */
#pragma GCC visibility push(default)

/* ===== Fortran helpers ======================================= */
void *hexpr_fortran_stdout(void);

/* ===== Versioning ============================================ */
const char *hexpr_version_string_(void);
void hexpr_version_major_(int *out);
void hexpr_version_minor_(int *out);
void hexpr_version_patch_(int *out);

/* ===== Status / error ======================================== */
const char *hexpr_last_error_(void);

/* ===== Lifecycle ============================================= */
void hexpr_init_(int *out_status);
void hexpr_shutdown_(void);

/* ===== DRT =================================================== */
hexpr_drt_t *hexpr_drt_create_(int *kind, int *nve, int *spin_x2,
                               int *ncore, int *nval, int *next);
hexpr_drt_t *hexpr_drt_from_csfset_(const int *norb,
                                    const unsigned char *csfs,
                                    const int *n_csfs);
void hexpr_drt_destroy_(hexpr_drt_t **drt);
void hexpr_drt_norb_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_nnodes_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_ncsfs_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_walk_ncsfs_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_kind_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_ncore_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_nval_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_next_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_nve_state_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_spin_x2_(const hexpr_drt_t **drt, int *out);
void hexpr_drt_csfs_(const hexpr_drt_t **drt, unsigned char *buf);
void hexpr_drt_walk_csfs_(const hexpr_drt_t **drt, unsigned char *buf);
void hexpr_drt_node_info_(const hexpr_drt_t **drt, int *node_index,
                          int *out_orb, int *out_nve, int *out_is,
                          int *out_status);
void hexpr_drt_show_(const hexpr_drt_t **drt, int *out_status);

/* ===== Expression ============================================ */
hexpr_expr_t *hexpr_expr_build_(const hexpr_drt_t **drt);
void hexpr_expr_destroy_(hexpr_expr_t **expr);
void hexpr_expr_nterms_full_(const hexpr_expr_t **expr, int64_t *out);
void hexpr_expr_nterms_one_(const hexpr_expr_t **expr, int64_t *out);
void hexpr_expr_read_(const hexpr_expr_t **expr, int *which,
                      int32_t *out_ij, int32_t *out_pqrs, double *out_coef,
                      int *out_status);
void hexpr_expr_show_full_(const hexpr_expr_t **expr, int *out_status);
void hexpr_expr_show_one_(const hexpr_expr_t **expr, int *out_status);

/* ===== CSF lookup ============================================ */
void hexpr_csf_index_(const hexpr_drt_t **drt, const signed char *step,
                      int *out_index, int *out_status);
void hexpr_csf_steps_(const hexpr_drt_t **drt, const int *index,
                      signed char *out_step, int *out_status);
void hexpr_step_parse_(const char *str, const int *norb,
                       signed char *out_step, int *out_status);

/* ===== CSF -> determinant expansion ========================== */
void hexpr_csf_ndets_(const hexpr_drt_t **drt, const int *index,
                      const int *ms_x2, int *out_ndets, int *out_status);
void hexpr_csf_to_dets_(const hexpr_drt_t **drt, const int *index,
                        const int *ms_x2, int *out_count,
                        uint64_t *out_alpha, uint64_t *out_beta,
                        double *out_coef, int *out_status);

/* ===== Per-pair / block / subspace =========================== */
void hexpr_pair_count_(const hexpr_drt_t **drt,
                       const int *bra, const int *ket,
                       const int *which, int *out_count, int *out_status);
void hexpr_pair_max_terms_(const hexpr_drt_t **drt,
                           const int *which, int *out_max, int *out_status);
void hexpr_pair_eval_(const hexpr_drt_t **drt,
                      const int *bra, const int *ket, const int *which,
                      int *out_count,
                      int32_t *out_ij, int32_t *out_pqrs, double *out_coef,
                      int *out_status);
void hexpr_block_count_(const hexpr_drt_t **drt,
                        const int *bra_indices, const int *n_bra,
                        const int *ket_indices, const int *n_ket,
                        const int *which, int *out_count, int *out_status);
void hexpr_block_max_terms_(const hexpr_drt_t **drt,
                            const int *n_bra, const int *n_ket,
                            const int *which, int *out_max, int *out_status);
void hexpr_block_eval_(const hexpr_drt_t **drt,
                       const int *bra_indices, const int *n_bra,
                       const int *ket_indices, const int *n_ket,
                       const int *which, int *out_count,
                       int32_t *out_bra_pos, int32_t *out_ket_pos,
                       int32_t *out_ij, int32_t *out_pqrs, double *out_coef,
                       int *out_status);
void hexpr_subspace_count_(const hexpr_drt_t **drt,
                           const int *indices, const int *n,
                           const int *which, int *out_count, int *out_status);
void hexpr_subspace_max_terms_(const hexpr_drt_t **drt,
                               const int *n, const int *which,
                               int *out_max, int *out_status);
void hexpr_subspace_eval_(const hexpr_drt_t **drt,
                          const int *indices, const int *n,
                          const int *which, int *out_count,
                          int32_t *out_bra_pos, int32_t *out_ket_pos,
                          int32_t *out_ij, int32_t *out_pqrs, double *out_coef,
                          int *out_status);

#pragma GCC visibility pop

#endif /* HEXPR_FORTRAN_H_ */
