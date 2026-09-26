/* csf/csf_pair.c -- public pair/block/subspace API via per-pair sub-DRT */
/**
 * @file csf_pair.c
 * @brief Public pair / block / subspace coupling API, evaluated one CSF pair at
 *        a time via a per-pair sub-DRT.
 *
 * This is the CSF-pairing half of the @c csf/ layer and a CONSUMER of the
 * @c csf_lookup.c canonical lookups and the @c drt.c accessors. The whole API
 * funnels through @ref pair_via_subdrt: given a (bra, ket) pair of canonical CSF
 * indices, it recovers their step paths (@ref hexpr_drt_walk_csf_steps), builds
 * a tiny sub-DRT containing just those one or two CSFs
 * (@ref hexpr_drt_from_csfset), evaluates the pair there
 * (@c pair_eval_internal in @c csf_pair_core.c), and relabels the local @c ij
 * back to the global pair index.
 *
 * INDEX SPACE. The public entry points take INPUT-ORDER CSF indices -- the
 * same space as @c hexpr_drt_csfs / @c hexpr_drt_ncsfs / @c hexpr_csf_steps --
 * and translate them to canonical (walk) indices through @ref input_to_walk
 * before anything else runs. The two spaces coincide for SOCI and FOCI and
 * differ only for a CSFSET sub-DRT, where the canonical walk may reorder the
 * input, may be a superset of it, or may dedup it.
 *
 * The emitted @c ij stays in the CANONICAL space: per @c docs/PQRS_INDEX_SPEC.md
 * (normative for v1.0.0) @c ij is @c triangle of the two Shavitt PATH WEIGHTS,
 * which is what the canonical index is. Input positions are not path weights,
 * so encoding them in @c ij would contradict that spec. Callers of the block
 * and subspace forms do not need @c ij to identify a pair: @c out_bra_pos and
 * @c out_ket_pos already give the position within the arrays they passed.
 *
 * The public surface is three families layered on that core:
 *   - pair: @ref hexpr_pair_count / @ref hexpr_pair_eval -- one (bra,ket) pair.
 *   - block: @ref hexpr_block_count / @ref hexpr_block_eval -- a bra-set x ket-set
 *     Cartesian product, both built on @ref block_eval_impl.
 *   - subspace: @ref hexpr_subspace_count / @ref hexpr_subspace_eval -- a block
 *     with identical bra and ket index sets.
 * Each family also has a @c *_max_terms upper-bound helper for buffer sizing.
 *
 * @note In-body comments describe only the data flow visible in the code. The
 *       GUGA addressing convention (canonical index <-> arc-weight sum) is OWNED
 *       by the @c csf_lookup.c todos and the one-/two-electron @c n_oel boundary
 *       by @c expr.c's @c hexpr_expr_build todo; those are NOT restated here.
 *       The only <tt>\@todo domain:</tt> originated here is the validity of the
 *       per-pair sub-DRT decomposition (see @ref pair_via_subdrt).
 */
#include "csf_pair_core.h"
#include "csf_lookup.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/** @brief Larger of two integers. @param a,b operands @return max(a,b). */
static int imax2(int a, int b) { return a > b ? a : b; }
/** @brief Smaller of two integers. @param a,b operands @return min(a,b). */
static int imin2(int a, int b) { return a < b ? a : b; }
/**
 * @brief Symmetric packed index of an unordered pair (a,b):
 *        @f$ \mathrm{max}(\mathrm{max}+1)/2 + \mathrm{min} @f$.
 * @param a,b the two CSF indices (order does not matter).
 * @return the packed pair index.
 *
 * Role: folds a (bra,ket) CSF-index pair into the single @c ij row address used
 * by the energy-expression terms. Here it converts both the sub-DRT-local pair
 * and the global pair to that address so @ref pair_via_subdrt can relabel one to
 * the other. (Same triangular packing convention used across the library; the
 * meaning of @c ij as a CSF-pair address is owned by @c expr.c's term layout.)
 */
static int triangle(int a, int b) { int mx=imax2(a,b),mn=imin2(a,b); return mx*(mx+1)/2+mn; }

/**
 * @brief Translate a caller-supplied INPUT-order CSF index to the canonical
 *        (walk) index the evaluator addresses.
 * @param drt the DRT.
 * @param idx input-order index, in @c [0, hexpr_drt_ncsfs()).
 * @param scratch caller-owned buffer of at least @c norb bytes, used to hold
 *        the row while it is looked up; its contents on return are unspecified.
 * @param[out] out_walk receives the canonical index.
 * @return @c HEXPR_OK, or @c HEXPR_ERR_RANGE if @p idx is out of range or its
 *         row has no canonical walk (which cannot happen for a DRT built from
 *         that row, and is reported rather than assumed away).
 *
 * Role: the single place where the public index space meets the internal one.
 * For SOCI and FOCI the two coincide and this is an identity in effect. For a
 * CSFSET sub-DRT the canonical walk may reorder, extend or dedup the input, so
 * the row is looked up by its step path rather than by position.
 *
 * @note Duplicate input rows (the dedup case) map to the SAME canonical index.
 *       That is correct: they are the same CSF, and the block forms still
 *       report their distinct positions through @c out_bra_pos / @c out_ket_pos.
 * @note The row is copied into @p scratch rather than passed straight out of
 *       the DRT-owned buffer, so this does not depend on whether the lookup
 *       takes its step argument as @c const.
 */
static hexpr_status_t input_to_walk(const hexpr_drt_t *drt, int idx,
                                    unsigned char *scratch, int *out_walk)
{
    if (idx < 0 || idx >= hexpr_drt_ncsfs(drt)) return HEXPR_ERR_RANGE;
    const unsigned char *rows = hexpr_drt_csfs(drt);
    if (!rows) return HEXPR_ERR_RANGE;
    int norb = hexpr_drt_norb(drt);
    memcpy(scratch, rows + (size_t)idx * (size_t)norb, (size_t)norb);
    if (hexpr_drt_walk_csf_index(drt, scratch, out_walk) != HEXPR_OK)
        return HEXPR_ERR_RANGE;
    return HEXPR_OK;
}

/* Evaluate one pair by building a per-pair sub-DRT, then patch global ij. */
/**
 * @brief Evaluate the coupling terms of one (bra, ket) CSF pair by building a
 *        sub-DRT that contains only those CSFs, then relabel @c ij to the global pair.
 * @param drt the full DRT the canonical indices @p bra / @p ket address.
 * @param bra,ket canonical CSF indices of the pair.
 * @return an ::ot_vec_t of contributing terms; the CALLER owns and must free
 *         @c .data. Empty (NULL data) on allocation failure or an invalid pair.
 *
 * Role: the single evaluation core every public pair/block/subspace call funnels
 * through. Data flow: (1) recover the two CSFs' step paths from the full DRT via
 * @ref hexpr_drt_walk_csf_steps; (2) pack them into a 1- or 2-row CSF buffer
 * (1 row when @c bra==ket) and build a per-pair sub-DRT with
 * @ref hexpr_drt_from_csfset; (3) re-index the same step paths INTO the sub-DRT
 * with @ref hexpr_drt_walk_csf_index to get their local indices; (4) evaluate
 * with @c pair_eval_internal; (5) destroy the sub-DRT and rewrite each term's
 * local pair address @c triangle(sub_bra,sub_ket) to the global
 * @c triangle(bra,ket).
 *
 * @note Every heap buffer (bra/ket steps, csf_buf, the sub-DRT) is released
 *       before return; only the returned @c result.data is handed to the caller.
 *       The indexing it performs follows the canonical convention owned by the
 *       @c csf_lookup.c todos and is not re-questioned here.
 *
 * DOMAIN (per paper, Honda & Noro): The pairwise reduction is exact because the
 * matrix element factorizes into per-orbital factors L_j (Eq.12/Eq.16, Sec. 2.4):
 * each L_j depends only on the bra step, ket step, and operator at orbital j (and
 * the adjacent spin labels), not on any other CSF in the space. The coupling of a
 * (bra,ket) pair is therefore a property of those two walks alone and is
 * independent of the ambient CSF set, so evaluating the pair inside an isolated
 * 1-2 CSF sub-DRT yields the same coefficients as within the full DRT. Building
 * that sub-DRT over the two supplied walks is the CSFSET construction of Sec. 4.1
 * (hexpr_drt_from_csfset). Canonical indexing follows the Shavitt addressing
 * documented in csf_lookup.c. See Sec. 2.4, Eq.12, Eq.16, Sec. 4.1.
 */
static ot_vec_t pair_via_subdrt(const hexpr_drt_t *drt, int bra, int ket)
{
    ot_vec_t empty = {NULL, 0, 0};
    int norb = hexpr_drt_norb(drt);
    /* Recover the two CSFs as step paths in the full DRT's canonical system. */
    unsigned char *bra_steps = malloc((size_t)norb);
    unsigned char *ket_steps = malloc((size_t)norb);
    if (!bra_steps || !ket_steps) { free(bra_steps); free(ket_steps); return empty; }
    /* A decode failure here must not fall through with a partly (or wholly)
     * uninitialized step buffer -- that buffer would otherwise go straight
     * into hexpr_drt_from_csfset below. */
    if (hexpr_drt_walk_csf_steps(drt, bra, bra_steps) != HEXPR_OK ||
        hexpr_drt_walk_csf_steps(drt, ket, ket_steps) != HEXPR_OK) {
        free(bra_steps); free(ket_steps);
        return empty;
    }

    /* Pack the pair into a tiny CSF set (one row if bra==ket) and build the
     * per-pair sub-DRT over it. */
    int n_sub = (bra == ket) ? 1 : 2;
    unsigned char *csf_buf = malloc((size_t)n_sub * (size_t)norb);
    if (!csf_buf) { free(bra_steps); free(ket_steps); return empty; }
    memcpy(csf_buf, bra_steps, (size_t)norb);
    if (n_sub == 2) memcpy(csf_buf + norb, ket_steps, (size_t)norb);

    hexpr_drt_t *sub = hexpr_drt_from_csfset(norb, csf_buf, n_sub);
    free(csf_buf);
    if (!sub) { free(bra_steps); free(ket_steps); return empty; }

    /* Re-index the same step paths within the sub-DRT (their local indices). */
    int sub_bra = 0, sub_ket = 0;
    if (hexpr_drt_walk_csf_index(sub, bra_steps, &sub_bra) != HEXPR_OK ||
        hexpr_drt_walk_csf_index(sub, ket_steps, &sub_ket) != HEXPR_OK) {
        free(bra_steps); free(ket_steps);
        hexpr_drt_destroy(sub);
        return empty;
    }
    free(bra_steps); free(ket_steps);

    ot_vec_t result = pair_eval_internal(sub, sub_bra, sub_ket);
    hexpr_drt_destroy(sub);

    /* Terms carry the sub-DRT-local pair address; rewrite it to the global one
     * so callers see ij in the full DRT's pair-index system. */
    int old_ij = triangle(sub_bra, sub_ket);
    int new_ij = triangle(bra, ket);
    for (int i = 0; i < result.size; i++)
        if (result.data[i].ij == old_ij) result.data[i].ij = new_ij;
    return result;
}

/**
 * @brief Count the coupling terms of one (bra, ket) pair, optionally
 *        one-electron only.
 * @param drt the DRT.
 * @param bra,ket INPUT-ORDER CSF indices, in @c [0, hexpr_drt_ncsfs());
 *        translated internally by @ref input_to_walk.
 * @param which @c HEXPR_EXPR_FULL counts all terms; @c HEXPR_EXPR_ONE counts only
 *              one-electron terms.
 * @param[out] out_count receives the term count.
 * @return @c HEXPR_OK; @c HEXPR_ERR_INVALID_ARG on NULL; @c HEXPR_ERR_RANGE if an
 *         index is out of @c [0, hexpr_drt_ncsfs()).
 *
 * Role: a thin count over @ref pair_via_subdrt. It evaluates the pair, then
 * tallies terms, skipping two-electron terms (@c pqrs >= n_oel) when @p which is
 * @c HEXPR_EXPR_ONE. The transient term vector is freed before return.
 *
 * @note Indices are INPUT-ORDER indices; see the INDEX SPACE note at the top of
 *       this file. The one-/two-electron split at
 *       @c n_oel = @f$ \mathit{norb}(\mathit{norb}+1)/2 @f$ is the same boundary
 *       documented by @c expr.c's @c hexpr_expr_build, and is not re-questioned here.
 */
hexpr_status_t hexpr_pair_count(const hexpr_drt_t *drt, int bra, int ket,
                                 hexpr_expr_which_t which, int *out_count)
{
    if (!drt || !out_count) return HEXPR_ERR_INVALID_ARG;
    unsigned char *scratch = malloc((size_t)hexpr_drt_norb(drt));
    if (!scratch) return HEXPR_ERR_OOM;
    int wbra = 0, wket = 0;
    if (input_to_walk(drt, bra, scratch, &wbra) != HEXPR_OK ||
        input_to_walk(drt, ket, scratch, &wket) != HEXPR_OK) {
        free(scratch); return HEXPR_ERR_RANGE;
    }
    free(scratch);

    int n_oel = hexpr_drt_norb(drt); n_oel = n_oel * (n_oel + 1) / 2;
    ot_vec_t v = pair_via_subdrt(drt, wbra, wket);
    /* Count terms, dropping two-electron ones (pqrs >= n_oel) in ONE mode. */
    int cnt = 0;
    for (int i = 0; i < v.size; i++) {
        if (which == HEXPR_EXPR_ONE && v.data[i].pqrs >= n_oel) continue;
        cnt++;
    }
    free(v.data);
    *out_count = cnt;
    return HEXPR_OK;
}

/**
 * @brief Evaluate one (bra, ket) pair, writing each term's @c ij / @c pqrs /
 *        @c coef into caller arrays.
 * @param drt the DRT.
 * @param bra,ket INPUT-ORDER CSF indices, in @c [0, hexpr_drt_ncsfs());
 *        translated internally by @ref input_to_walk.
 * @param which @c HEXPR_EXPR_FULL emits all terms; @c HEXPR_EXPR_ONE emits only
 *              one-electron terms.
 * @param[out] out_count receives the number of terms written.
 * @param[out] out_ij,out_pqrs,out_coef parallel output arrays; each may be NULL
 *             to skip that column.
 * @return @c HEXPR_OK; @c HEXPR_ERR_INVALID_ARG on NULL drt/out_count;
 *         @c HEXPR_ERR_RANGE if an index is out of range.
 *
 * Role: the emitting twin of @ref hexpr_pair_count -- same evaluate-and-filter
 * over @ref pair_via_subdrt, but it copies the surviving terms into the parallel
 * output arrays (writing column @c k for each kept term) instead of only
 * counting. Size the arrays with @ref hexpr_pair_max_terms.
 *
 * @note Input-order indices; the emitted @c ij remains a CANONICAL pair address
 *       (see the INDEX SPACE note at the top of this file). @c n_oel
 *       one-/two-electron split per @c expr.c (not re-questioned here). Each
 *       output pointer is optional.
 */
hexpr_status_t hexpr_pair_eval(const hexpr_drt_t *drt, int bra, int ket,
                                hexpr_expr_which_t which, int *out_count,
                                int32_t *out_ij, int32_t *out_pqrs, double *out_coef)
{
    if (!drt || !out_count) return HEXPR_ERR_INVALID_ARG;
    unsigned char *scratch = malloc((size_t)hexpr_drt_norb(drt));
    if (!scratch) return HEXPR_ERR_OOM;
    int wbra = 0, wket = 0;
    if (input_to_walk(drt, bra, scratch, &wbra) != HEXPR_OK ||
        input_to_walk(drt, ket, scratch, &wket) != HEXPR_OK) {
        free(scratch); return HEXPR_ERR_RANGE;
    }
    free(scratch);

    int n_oel = hexpr_drt_norb(drt); n_oel = n_oel * (n_oel + 1) / 2;
    ot_vec_t v = pair_via_subdrt(drt, wbra, wket);
    /* Copy kept terms into the (optional) parallel output columns; k counts
     * survivors, skipping two-electron terms in ONE mode. */
    int k = 0;
    for (int i = 0; i < v.size; i++) {
        const out_term_t *t = &v.data[i];
        if (which == HEXPR_EXPR_ONE && t->pqrs >= n_oel) continue;
        if (out_ij)   out_ij[k]   = (int32_t)t->ij;
        if (out_pqrs) out_pqrs[k] = (int32_t)t->pqrs;
        if (out_coef) out_coef[k] = t->coef;
        k++;
    }
    free(v.data);
    *out_count = k;
    return HEXPR_OK;
}

/**
 * @brief Shared engine for block count/eval: evaluate every (bra,ket) pair in a
 *        bra-set x ket-set Cartesian product, counting or emitting terms.
 * @param drt the DRT.
 * @param bra_indices,n_bra the bra INPUT-ORDER index set and its size.
 * @param ket_indices,n_ket the ket INPUT-ORDER index set and its size.
 * @param which @c HEXPR_EXPR_FULL / @c HEXPR_EXPR_ONE term filter.
 * @param[out] out_count receives the total number of terms (counted or emitted).
 * @param[out] out_bra_pos,out_ket_pos per-term SET POSITIONS @c ip / @c jp (not
 *             the CSF indices); each optional.
 * @param[out] out_ij,out_pqrs,out_coef per-term columns; each optional.
 * @param emit 0 = count only (no writes); 1 = also write the output columns.
 * @return @c HEXPR_OK; @c HEXPR_ERR_INVALID_ARG on NULL/negative sizes/missing
 *         index arrays; @c HEXPR_ERR_OOM if the translation buffers cannot be
 *         allocated; @c HEXPR_ERR_RANGE if any index is out of range.
 *
 * Role: the single implementation behind @ref hexpr_block_count (@c emit=0) and
 * @ref hexpr_block_eval (@c emit=1); the subspace calls reach it through those
 * with identical bra/ket sets. It translates both index sets to canonical
 * indices up front -- which validates them, so the evaluation loop cannot fault
 * -- then runs a nested @c ip x @c jp loop, evaluates each pair with
 * @ref pair_via_subdrt, and walks the resulting terms with the same
 * one-electron filter as the pair calls.
 * The running counter @c k is the write cursor when emitting and simply the
 * total when counting, so both modes share one code path.
 *
 * @note The emitted @c out_bra_pos / @c out_ket_pos are positions WITHIN the
 *       arrays the caller passed, letting a caller place each term back into a
 *       dense block matrix. They, not @c ij, are how a term is attributed to a
 *       pair here; @c ij stays a canonical pair address (see the INDEX SPACE
 *       note at the top of this file). Per-pair correctness and the @c n_oel
 *       split are owned upstream (@ref pair_via_subdrt and @c expr.c
 *       respectively) and not restated here.
 */
static hexpr_status_t block_eval_impl(
    const hexpr_drt_t *drt,
    const int *bra_indices, int n_bra,
    const int *ket_indices, int n_ket,
    hexpr_expr_which_t which, int *out_count,
    int32_t *out_bra_pos, int32_t *out_ket_pos,
    int32_t *out_ij, int32_t *out_pqrs, double *out_coef,
    int emit)
{
    if (!drt || !out_count) return HEXPR_ERR_INVALID_ARG;
    if (n_bra < 0 || n_ket < 0) return HEXPR_ERR_INVALID_ARG;
    if ((n_bra > 0 && !bra_indices) || (n_ket > 0 && !ket_indices))
        return HEXPR_ERR_INVALID_ARG;

    /* Translate both sets to canonical indices up front. This also validates
     * them, so the evaluation loop below cannot fault on a bad index. */
    int *bra_walk = NULL, *ket_walk = NULL;
    unsigned char *scratch = malloc((size_t)hexpr_drt_norb(drt));
    if (!scratch) return HEXPR_ERR_OOM;
    if (n_bra > 0) {
        bra_walk = malloc((size_t)n_bra * sizeof(int));
        if (!bra_walk) { free(scratch); return HEXPR_ERR_OOM; }
    }
    if (n_ket > 0) {
        ket_walk = malloc((size_t)n_ket * sizeof(int));
        if (!ket_walk) { free(scratch); free(bra_walk); return HEXPR_ERR_OOM; }
    }
    for (int ip = 0; ip < n_bra; ip++)
        if (input_to_walk(drt, bra_indices[ip], scratch, &bra_walk[ip]) != HEXPR_OK) {
            free(scratch); free(bra_walk); free(ket_walk); return HEXPR_ERR_RANGE;
        }
    for (int jp = 0; jp < n_ket; jp++)
        if (input_to_walk(drt, ket_indices[jp], scratch, &ket_walk[jp]) != HEXPR_OK) {
            free(scratch); free(bra_walk); free(ket_walk); return HEXPR_ERR_RANGE;
        }
    free(scratch);

    int n_oel = hexpr_drt_norb(drt); n_oel = n_oel * (n_oel + 1) / 2;
    /* k is the shared cursor: a write index when emit=1, else just a tally. */
    int k = 0;
    for (int ip = 0; ip < n_bra; ip++) {
        for (int jp = 0; jp < n_ket; jp++) {
            ot_vec_t v = pair_via_subdrt(drt, bra_walk[ip], ket_walk[jp]);
            for (int i = 0; i < v.size; i++) {
                const out_term_t *t = &v.data[i];
                if (which == HEXPR_EXPR_ONE && t->pqrs >= n_oel) continue;
                /* Emit mode also records the set positions ip/jp alongside the
                 * term so the caller can locate it within the block. */
                if (emit) {
                    if (out_bra_pos) out_bra_pos[k] = (int32_t)ip;
                    if (out_ket_pos) out_ket_pos[k] = (int32_t)jp;
                    if (out_ij)      out_ij[k]      = (int32_t)t->ij;
                    if (out_pqrs)    out_pqrs[k]    = (int32_t)t->pqrs;
                    if (out_coef)    out_coef[k]    = t->coef;
                }
                k++;
            }
            free(v.data);
        }
    }
    free(bra_walk); free(ket_walk);
    *out_count = k;
    return HEXPR_OK;
}

/**
 * @brief Count all coupling terms over a bra-set x ket-set block.
 * @param drt the DRT.
 * @param bra_indices,n_bra bra INPUT-ORDER index set and size.
 * @param ket_indices,n_ket ket INPUT-ORDER index set and size.
 * @param which term filter (@c HEXPR_EXPR_FULL / @c HEXPR_EXPR_ONE).
 * @param[out] out_count receives the total term count.
 * @return status from @ref block_eval_impl.
 * @note Thin count-only wrapper (@c emit=0) over @ref block_eval_impl; use its
 *       result to size the arrays for @ref hexpr_block_eval.
 */
hexpr_status_t hexpr_block_count(
    const hexpr_drt_t *drt,
    const int *bra_indices, int n_bra,
    const int *ket_indices, int n_ket,
    hexpr_expr_which_t which, int *out_count)
{
    return block_eval_impl(drt, bra_indices, n_bra, ket_indices, n_ket,
                            which, out_count, NULL, NULL, NULL, NULL, NULL, 0);
}

/**
 * @brief Evaluate all coupling terms over a bra-set x ket-set block.
 * @param drt the DRT.
 * @param bra_indices,n_bra bra INPUT-ORDER index set and size.
 * @param ket_indices,n_ket ket INPUT-ORDER index set and size.
 * @param which term filter (@c HEXPR_EXPR_FULL / @c HEXPR_EXPR_ONE).
 * @param[out] out_count receives the number of terms written.
 * @param[out] out_bra_pos,out_ket_pos per-term set positions (each optional).
 * @param[out] out_ij,out_pqrs,out_coef per-term columns (each optional).
 * @return status from @ref block_eval_impl.
 * @note Thin emit wrapper (@c emit=1) over @ref block_eval_impl; size outputs
 *       with @ref hexpr_block_max_terms or an exact @ref hexpr_block_count.
 */
hexpr_status_t hexpr_block_eval(
    const hexpr_drt_t *drt,
    const int *bra_indices, int n_bra,
    const int *ket_indices, int n_ket,
    hexpr_expr_which_t which, int *out_count,
    int32_t *out_bra_pos, int32_t *out_ket_pos,
    int32_t *out_ij, int32_t *out_pqrs, double *out_coef)
{
    return block_eval_impl(drt, bra_indices, n_bra, ket_indices, n_ket,
                            which, out_count,
                            out_bra_pos, out_ket_pos, out_ij, out_pqrs, out_coef, 1);
}

/**
 * @brief Count coupling terms over a subspace (a block with one index set used
 *        as both bra and ket).
 * @param drt the DRT.
 * @param indices,n the INPUT-ORDER index set and its size.
 * @param which term filter.
 * @param[out] out_count receives the total term count.
 * @return status from @ref hexpr_block_count.
 * @note Delegates to @ref hexpr_block_count with @p indices as both bra and ket.
 */
hexpr_status_t hexpr_subspace_count(
    const hexpr_drt_t *drt, const int *indices, int n,
    hexpr_expr_which_t which, int *out_count)
{
    return hexpr_block_count(drt, indices, n, indices, n, which, out_count);
}

/**
 * @brief Evaluate coupling terms over a subspace (one index set as both bra and ket).
 * @param drt the DRT.
 * @param indices,n the INPUT-ORDER index set and its size.
 * @param which term filter.
 * @param[out] out_count receives the number of terms written.
 * @param[out] out_bra_pos,out_ket_pos per-term set positions (each optional).
 * @param[out] out_ij,out_pqrs,out_coef per-term columns (each optional).
 * @return status from @ref hexpr_block_eval.
 * @note Delegates to @ref hexpr_block_eval with @p indices as both bra and ket.
 */
hexpr_status_t hexpr_subspace_eval(
    const hexpr_drt_t *drt, const int *indices, int n,
    hexpr_expr_which_t which, int *out_count,
    int32_t *out_bra_pos, int32_t *out_ket_pos,
    int32_t *out_ij, int32_t *out_pqrs, double *out_coef)
{
    return hexpr_block_eval(drt, indices, n, indices, n, which, out_count,
                             out_bra_pos, out_ket_pos, out_ij, out_pqrs, out_coef);
}

/* ===== Max-terms upper bounds ============================== */

/**
 * @brief Upper bound on coupling terms a single pair can produce.
 * @param norb orbital count.
 * @param which @c HEXPR_EXPR_ONE bounds one-electron terms only; otherwise both.
 * @return @c n_oel for ONE, else @c n_oel + n_oel*(n_oel+1)/2, with
 *         @c n_oel = @f$ \mathit{norb}(\mathit{norb}+1)/2 @f$ (one-electron count
 *         plus the count of two-electron index pairs).
 * @note A worst-case sizing bound, not an exact count; shared by the three
 *       public @c *_max_terms helpers.
 */
static int per_pair_max_terms(int norb, hexpr_expr_which_t which)
{
    int n_oel = norb * (norb + 1) / 2;
    if (which == HEXPR_EXPR_ONE) return n_oel;
    return n_oel + n_oel * (n_oel + 1) / 2;
}

/**
 * @brief Upper bound on terms from a single pair, for buffer sizing.
 * @param drt the DRT.
 * @param which term filter.
 * @param[out] out_max receives the bound.
 * @return @c HEXPR_OK, or @c HEXPR_ERR_INVALID_ARG on NULL.
 * @note Wraps @ref per_pair_max_terms at this DRT's @c norb; size
 *       @ref hexpr_pair_eval outputs with it.
 */
hexpr_status_t hexpr_pair_max_terms(const hexpr_drt_t *drt,
                                     hexpr_expr_which_t which, int *out_max)
{
    if (!drt || !out_max) return HEXPR_ERR_INVALID_ARG;
    *out_max = per_pair_max_terms(hexpr_drt_norb(drt), which);
    return HEXPR_OK;
}

/**
 * @brief Upper bound on terms from a bra-set x ket-set block, for buffer sizing.
 * @param drt the DRT.
 * @param n_bra,n_ket block dimensions.
 * @param which term filter.
 * @param[out] out_max receives the bound.
 * @return @c HEXPR_OK; @c HEXPR_ERR_INVALID_ARG on NULL or negative sizes;
 *         @c HEXPR_ERR_RANGE if @c n_bra*n_ket*per_pair overflows @c int.
 * @note Computes the product in @c int64_t and rejects overflow before casting
 *       back to @c int (per_pair from @ref per_pair_max_terms).
 */
hexpr_status_t hexpr_block_max_terms(const hexpr_drt_t *drt,
                                      int n_bra, int n_ket,
                                      hexpr_expr_which_t which, int *out_max)
{
    if (!drt || !out_max) return HEXPR_ERR_INVALID_ARG;
    if (n_bra < 0 || n_ket < 0) return HEXPR_ERR_INVALID_ARG;
    int per_pair = per_pair_max_terms(hexpr_drt_norb(drt), which);
    int64_t product = (int64_t)n_bra * (int64_t)n_ket * (int64_t)per_pair;
    if (product > (int64_t)INT_MAX) return HEXPR_ERR_RANGE;
    *out_max = (int)product;
    return HEXPR_OK;
}

/**
 * @brief Upper bound on terms from an n x n subspace, for buffer sizing.
 * @param drt the DRT.
 * @param n subspace dimension (used as both bra and ket size).
 * @param which term filter.
 * @param[out] out_max receives the bound.
 * @return @c HEXPR_OK; @c HEXPR_ERR_INVALID_ARG on NULL or negative @p n;
 *         @c HEXPR_ERR_RANGE if @c n*n*per_pair overflows @c int.
 * @note The @ref hexpr_block_max_terms case with @c n_bra==n_ket==n; same
 *       @c int64_t overflow guard.
 */
hexpr_status_t hexpr_subspace_max_terms(const hexpr_drt_t *drt,
                                         int n,
                                         hexpr_expr_which_t which, int *out_max)
{
    if (!drt || !out_max) return HEXPR_ERR_INVALID_ARG;
    if (n < 0) return HEXPR_ERR_INVALID_ARG;
    int per_pair = per_pair_max_terms(hexpr_drt_norb(drt), which);
    int64_t product = (int64_t)n * (int64_t)n * (int64_t)per_pair;
    if (product > (int64_t)INT_MAX) return HEXPR_ERR_RANGE;
    *out_max = (int)product;
    return HEXPR_OK;
}
