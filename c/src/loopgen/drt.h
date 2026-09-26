/* loopgen/drt.h -- Distinct Row Table (DRT) internal header */
#ifndef HEXPR_DRT_H_
#define HEXPR_DRT_H_

#include <stdio.h>
#include <stdint.h>

/* hexpr_drt_t, hexpr_drt_kind_t, the two builders and the scalar accessors
 * are public API and come from the installed header.
 *
 * They used to be duplicated below, behind #ifndef HEXPR_H_, so that a
 * translation unit including only this file still saw them. drt.c is such a
 * translation unit, which made that copy -- not hexpr.h -- the declaration
 * the compiler saw when compiling their definitions. Two consequences, both
 * observed: every public function defined here was hidden when the library
 * was first built with -fvisibility=hidden, because the pragma was in
 * hexpr.h and the copy had none; and the copy had drifted, declaring
 * hexpr_drt_show and hexpr_drt_node_info as returning int rather than
 * hexpr_status_t, and documenting hexpr_drt_show's error as HEXPR_ERR_IO
 * when the code returns HEXPR_ERR_INVALID_ARG.
 *
 * One declaration cannot drift from itself. Do not reintroduce a copy: add
 * public declarations to hexpr.h only. */
#include "hexpr.h"

/* ------------------------------------------------------------------ */
/* Internal only: not declared in hexpr.h, hidden in libhexpr.so       */
/* ------------------------------------------------------------------ */

/* hexpr_drt_free is internal (hexpr_drt_destroy is the public alias) */
void hexpr_drt_free(hexpr_drt_t *drt);

/* Per-node accessors */
int hexpr_drt_orb  (const hexpr_drt_t *drt, int inode);
int hexpr_drt_nve  (const hexpr_drt_t *drt, int inode);
int hexpr_drt_spin (const hexpr_drt_t *drt, int inode);   /* (2S+1) */
int hexpr_drt_upper_arc (const hexpr_drt_t *drt, int inode, int w);
int hexpr_drt_lower_arc (const hexpr_drt_t *drt, int inode, int w);
int hexpr_drt_arc_weight(const hexpr_drt_t *drt, int inode, int w);
int hexpr_drt_number_upper_arcs(const hexpr_drt_t *drt, int inode);
int hexpr_drt_number_lower_arcs(const hexpr_drt_t *drt, int inode);

/* Path weight arrays -- do NOT free; valid until hexpr_drt_free(). */
const int *hexpr_drt_lower_path_weights(const hexpr_drt_t *drt, int inode, int *out_n);
const int *hexpr_drt_upper_path_weights(const hexpr_drt_t *drt, int inode, int *out_n);

/* CSF list -- caller frees each element. Returns ncsf, or -1 on error. */
int hexpr_drt_mk_csf_list(const hexpr_drt_t *drt, char ***out_steps);

#endif /* HEXPR_DRT_H_ */
