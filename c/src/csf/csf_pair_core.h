/* csf/csf_pair_core.h -- internal per-pair Brooks core */
#ifndef HEXPR_CSF_PAIR_CORE_H_
#define HEXPR_CSF_PAIR_CORE_H_

#include "hexpr.h"
#include "expr_internal.h"

/* Core per-pair evaluation. drt may be a sub-DRT or a full DRT.
 * (sub_bra, sub_ket) are path-weight indices within drt's CSF list.
 * Returns all contributing terms; caller must free result.data. */
ot_vec_t pair_eval_internal(const hexpr_drt_t *drt, int sub_bra, int sub_ket);

#endif /* HEXPR_CSF_PAIR_CORE_H_ */
