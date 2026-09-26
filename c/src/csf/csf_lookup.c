/* csf/csf_lookup.c -- CSF index/step-code lookup implementation */
/**
 * @file csf_lookup.c
 * @brief CSF locating/indexing: map between a step-code path and its CSF index,
 *        in both the canonical GUGA-walk system and the user-facing input system.
 *
 * Given a DRT and its canonical CSF walks (built in @c drt.c), this unit answers
 * "what is the index of this CSF?" and "what is the CSF at this index?". It has
 * two layers:
 *   - The internal CANONICAL lookups @ref hexpr_drt_walk_csf_index and
 *     @ref hexpr_drt_walk_csf_steps address the canonical (possibly superset)
 *     walk list for every DRT kind by walking the graph through the
 *     @c hexpr_drt_* accessors (upper arcs, arc weights, upper-arc counts).
 *     These are the load-bearing lookups that @c csf_pair.c relies on.
 *   - The public USER-FACING lookups @ref hexpr_csf_index / @ref hexpr_csf_steps
 *     speak the retained input set for a CSFSET sub-DRT (input order), and
 *     delegate to the canonical lookups for SOCI/FOCI (where no input set exists).
 *   - @ref hexpr_step_parse converts an "eudf" string into the step-code bytes
 *     the lookups consume.
 */
#include "csf_lookup.h"
#include "drt.h"
#include <string.h>

/* ================================================================== */
/* Internal CANONICAL lookups (GUGA walk order).  NOT in the public    */
/* header.  These are the load-bearing lookups the evaluation path     */
/* (csf_pair.c) relies on, and they address the canonical (possibly    */
/* superset) CSF list for ALL kinds including CSFSET sub-DRTs.          */
/* ================================================================== */

/* Walk top-down from root (node 0) following step[k] at each orbital level.
 * CSF index = sum of arc_weight[step[k]] along the path (canonical index). */
/**
 * @brief Map a step-code path to its canonical CSF index by walking the DRT.
 * @param drt       the DRT to walk.
 * @param step      array of @c norb step codes (each 0..3, one per orbital level).
 * @param[out] out_index receives the canonical index (untouched on error).
 * @return @c HEXPR_OK, or @c HEXPR_ERR_INVALID_ARG for a NULL argument, a step
 *         byte > 3, or a step that has no arc at the current node.
 *
 * Role: the forward canonical lookup, the inverse of
 * @ref hexpr_drt_walk_csf_steps and the addressing primitive @c csf_pair.c uses.
 * Starting at the root (node 0), it follows @c upper_arc for @c step[k] at each
 * level via @ref hexpr_drt_upper_arc, accumulating @ref hexpr_drt_arc_weight
 * into a running @c index; the index is the sum of arc weights along the path.
 * A missing arc (@c child < 0) means the path is not a valid walk.
 *
 * DOMAIN (Shavitt DRT; the paper uses this addressing but does not re-derive it):
 * The canonical index is the sum of arc weights because the arc weights encode
 * a reverse-lexical ordering of walks (Shavitt, Int.J.Quant.Chem. S11/S12,
 * refs [4,5]). Each arc weight counts the walks that, at this level, branch to a
 * lexically-earlier arc; summing them along a path yields the number of walks
 * ordered before it, i.e. its position (canonical index) in the walk list. The
 * paper builds the DRT "following Shavitt" (Sec. 4.1) but does not restate the
 * addressing formula, so the convention is cited to Shavitt rather than to an
 * equation here. Set up by store_number_upper_arcs in drt.c.
 */
hexpr_status_t hexpr_drt_walk_csf_index(const hexpr_drt_t *drt,
                                         const unsigned char *step,
                                         int *out_index)
{
    if (!drt || !step || !out_index) return HEXPR_ERR_INVALID_ARG;
    int norb = hexpr_drt_norb(drt);
    int node = 0;
    int index = 0;
    /* Descend one orbital level per step, accumulating arc weights into index
     * and advancing to the arc's child node. */
    for (int k = 0; k < norb; k++) {
        unsigned char s = step[k];
        if (s > 3) return HEXPR_ERR_INVALID_ARG;
        int child = hexpr_drt_upper_arc(drt, node, (int)s);
        if (child < 0) return HEXPR_ERR_INVALID_ARG;
        index += hexpr_drt_arc_weight(drt, node, (int)s);
        node = child;
    }
    *out_index = index;
    return HEXPR_OK;
}

/* Inverse canonical walk: given a CANONICAL index, recover the step-code path.
 * At each node, find the arc w such that:
 *   arc_weight[w] <= remaining < arc_weight[w] + max(nua[child], 1)
 * Then subtract arc_weight[w] from remaining and move to child. */
/**
 * @brief Recover the step-code path of a canonical CSF index by walking the DRT.
 * @param drt       the DRT to walk.
 * @param index     canonical index in @c [0, walk_ncsfs).
 * @param[out] out_step receives @c norb step codes (0..3), one per level.
 * @return @c HEXPR_OK; @c HEXPR_ERR_INVALID_ARG on NULL; @c HEXPR_ERR_RANGE if
 *         @p index is out of range; @c HEXPR_ERR_INTERNAL if no arc matches
 *         (which would indicate a malformed DRT).
 *
 * Role: the inverse of @ref hexpr_drt_walk_csf_index. From the root it greedily
 * picks, at each level, the arc whose weight interval contains the remaining
 * index -- i.e. the arc @c w with
 * @f$ \mathit{aw}_w \le \mathit{remaining} < \mathit{aw}_w + \max(\mathit{nua}_{child},1) @f$,
 * where @c nua is @ref hexpr_drt_number_upper_arcs (the child's CSF count, the
 * "width" of its slot). It then subtracts the arc weight and descends. The
 * @c max(nua,1) guard accounts for the top node, whose @c nua is 0 yet still
 * represents one CSF (per the pre-existing comment).
 *
 * DOMAIN (Shavitt DRT; addressing used by the paper, not re-derived there):
 * Each arc owns the index sub-range [arc_weight, arc_weight + subtree_csf_count)
 * because the reverse-lexical ordering places all walks through a given arc in a
 * contiguous block: arc_weight is the block's offset (walks through earlier arcs)
 * and the child's CSF count (number_upper_arcs) is its width. The per-arc blocks
 * therefore tile the canonical index space [0, walk_ncsfs) exactly, which is what
 * makes the greedy inverse walk here well-defined. This is the Shavitt reverse-
 * lexical addressing (refs [4,5]); the paper references the DRT construction in
 * Sec. 4.1 without restating the interval property. The max(nua,1) guard handles
 * the top node (nua 0 but one CSF).
 */
hexpr_status_t hexpr_drt_walk_csf_steps(const hexpr_drt_t *drt,
                                         int index,
                                         unsigned char *out_step)
{
    if (!drt || !out_step) return HEXPR_ERR_INVALID_ARG;
    int ncsfs = hexpr_drt_walk_ncsfs(drt);   /* canonical count */
    if (index < 0 || index >= ncsfs) return HEXPR_ERR_RANGE;
    int norb = hexpr_drt_norb(drt);
    int node = 0;
    int remaining = index;
    /* At each level, scan the four arcs and take the one whose [aw, aw+sn) slot
     * contains `remaining`; subtract its base weight and descend into the child. */
    for (int k = 0; k < norb; k++) {
        int taken = 0;
        for (int w = 0; w < 4; w++) {
            int child = hexpr_drt_upper_arc(drt, node, w);
            if (child < 0) continue;
            int aw = hexpr_drt_arc_weight(drt, node, w);
            int sn = hexpr_drt_number_upper_arcs(drt, child);
            if (sn < 1) sn = 1;  /* top node: nua==0 but represents 1 CSF */
            if (remaining < aw + sn) {
                out_step[k] = (unsigned char)w;
                remaining -= aw;
                node = child;
                taken = 1;
                break;
            }
        }
        if (!taken) return HEXPR_ERR_INTERNAL;
    }
    return HEXPR_OK;
}

/* ================================================================== */
/* Public USER-FACING lookups (input system for CSFSET sub-DRTs).      */
/* For SOCI/FOCI there is no input set, so they delegate to the        */
/* canonical lookups (canonical == the user's parametric result).      */
/* ================================================================== */

/* For a CSFSET sub-DRT: return the 0-based position of `step` within the
 * retained input set (input order). A `step` sequence that is valid in the
 * DRT but is NOT one of the N input rows (e.g. a superset-only canonical
 * CSF) has no input-system index and returns HEXPR_ERR_RANGE. Malformed
 * step bytes (>3) return HEXPR_ERR_INVALID_ARG. For SOCI/FOCI: the
 * canonical index (unchanged behavior). */
/**
 * @brief User-facing CSF index: position in the input set for CSFSET, else the
 *        canonical index.
 * @param drt       the DRT.
 * @param step      @c norb step codes (0..3) identifying the CSF.
 * @param[out] out_index receives the index (untouched on error).
 * @return @c HEXPR_OK; @c HEXPR_ERR_INVALID_ARG on NULL or a step byte > 3;
 *         @c HEXPR_ERR_RANGE for a CSFSET when @p step is a valid walk but not
 *         one of the retained input rows.
 *
 * Role: the public counterpart to the canonical
 * @ref hexpr_drt_walk_csf_index lookup.
 * For a CSFSET sub-DRT it validates the step bytes then LINEARLY SEARCHES the
 * retained input rows (@ref hexpr_drt_csfs) by @c memcmp, returning the 0-based
 * input position -- so the result is purely the caller's input ordering, with no
 * GUGA convention involved. For SOCI/FOCI it delegates unchanged to
 * @ref hexpr_drt_walk_csf_index.
 *
 * @note Borrows the input rows from @p drt (no copy); valid until the DRT is
 *       freed. O(N*norb) for the CSFSET search.
 */
hexpr_status_t hexpr_csf_index(const hexpr_drt_t *drt,
                                const unsigned char *step,
                                int *out_index)
{
    if (!drt || !step || !out_index) return HEXPR_ERR_INVALID_ARG;
    /* CSFSET: the user system is the input order, so search the input rows;
     * any other kind has no input set and falls through to the canonical walk. */
    if (hexpr_drt_kind(drt) == HEXPR_DRT_CSFSET) {
        int norb = hexpr_drt_norb(drt);
        for (int k = 0; k < norb; k++)
            if (step[k] > 3) return HEXPR_ERR_INVALID_ARG;
        int n = hexpr_drt_ncsfs(drt);                 /* input count N */
        const unsigned char *rows = hexpr_drt_csfs(drt);  /* input rows */
        /* Linear scan: first row equal to `step` gives its input position. */
        for (int i = 0; i < n; i++) {
            if (memcmp(rows + (size_t)i * (size_t)norb, step,
                       (size_t)norb) == 0) {
                *out_index = i;
                return HEXPR_OK;
            }
        }
        return HEXPR_ERR_RANGE;   /* valid steps, but not an input CSF */
    }
    return hexpr_drt_walk_csf_index(drt, step, out_index);
}

/* For a CSFSET sub-DRT: return input row `index` (0..N-1, input order). For
 * SOCI/FOCI: the canonical step path for the canonical index (unchanged). */
/**
 * @brief User-facing inverse lookup: input row @p index for CSFSET, else the
 *        canonical step path.
 * @param drt       the DRT.
 * @param index     index in @c [0, N) for CSFSET (input order), else canonical.
 * @param[out] out_step receives @c norb step codes (0..3).
 * @return @c HEXPR_OK; @c HEXPR_ERR_INVALID_ARG on NULL; @c HEXPR_ERR_RANGE if
 *         @p index is out of range.
 *
 * Role: the inverse of @ref hexpr_csf_index and the public counterpart to
 * @ref hexpr_drt_walk_csf_steps. For a CSFSET sub-DRT it simply @c memcpy's the
 * requested input row out of @ref hexpr_drt_csfs (a direct row copy in input
 * order, no walk). For SOCI/FOCI it delegates unchanged to
 * @ref hexpr_drt_walk_csf_steps.
 *
 * @note @p out_step must have room for @c norb bytes.
 */
hexpr_status_t hexpr_csf_steps(const hexpr_drt_t *drt,
                                int index,
                                unsigned char *out_step)
{
    if (!drt || !out_step) return HEXPR_ERR_INVALID_ARG;
    /* CSFSET: copy input row `index` directly; otherwise walk canonically. */
    if (hexpr_drt_kind(drt) == HEXPR_DRT_CSFSET) {
        int n = hexpr_drt_ncsfs(drt);                 /* input count N */
        if (index < 0 || index >= n) return HEXPR_ERR_RANGE;
        int norb = hexpr_drt_norb(drt);
        memcpy(out_step, hexpr_drt_csfs(drt) + (size_t)index * (size_t)norb,
               (size_t)norb);
        return HEXPR_OK;
    }
    return hexpr_drt_walk_csf_steps(drt, index, out_step);
}

/* Parse "eudf..." string (case-insensitive) into step-code byte array. */
/**
 * @brief Parse an "eudf" step string (case-insensitive) into step-code bytes.
 * @param str       NUL-terminated string of at least @p norb step characters
 *                  from {e,u,d,f} (any case).
 * @param norb      number of characters to parse.
 * @param[out] out_step receives @p norb step codes (@c HEXPR_STEP_E/U/D/F).
 * @return @c HEXPR_OK, or @c HEXPR_ERR_INVALID_ARG on NULL, @c norb<=0, a short
 *         string (early NUL), or any character outside {e,u,d,f}.
 * @note Produces the step-code array consumed by @ref hexpr_csf_index /
 *       @ref hexpr_drt_walk_csf_index. Reads exactly @p norb characters.
 */
hexpr_status_t hexpr_step_parse(const char *str, int norb,
                                 unsigned char *out_step)
{
    if (!str || !out_step || norb <= 0) return HEXPR_ERR_INVALID_ARG;
    for (int k = 0; k < norb; k++) {
        char c = str[k];
        if (c == '\0') return HEXPR_ERR_INVALID_ARG;
        switch (c) {
            case 'e': case 'E': out_step[k] = HEXPR_STEP_E; break;
            case 'u': case 'U': out_step[k] = HEXPR_STEP_U; break;
            case 'd': case 'D': out_step[k] = HEXPR_STEP_D; break;
            case 'f': case 'F': out_step[k] = HEXPR_STEP_F; break;
            default: return HEXPR_ERR_INVALID_ARG;
        }
    }
    return HEXPR_OK;
}
