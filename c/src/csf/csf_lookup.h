/* csf/csf_lookup.h -- CSF index/step-code lookup API */
#ifndef HEXPR_CSF_LOOKUP_H_
#define HEXPR_CSF_LOOKUP_H_

#include "hexpr.h"

/* User-facing lookups (hexpr_csf_index / hexpr_csf_steps / hexpr_step_parse)
 * are declared in the public hexpr.h, included above. */

/* Internal CANONICAL (GUGA-walk) lookups -- NOT exported in the public
 * header. The evaluation path (csf_pair.c) relies on these to address the
 * canonical (possibly superset) CSF list for any kind, including CSFSET. */
hexpr_status_t hexpr_drt_walk_csf_index(const hexpr_drt_t *drt,
                                         const unsigned char *step,
                                         int *out_index);

hexpr_status_t hexpr_drt_walk_csf_steps(const hexpr_drt_t *drt,
                                         int index,
                                         unsigned char *out_step);

#endif /* HEXPR_CSF_LOOKUP_H_ */
