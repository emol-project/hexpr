/* csf/csf_pair_core.c -- per-pair Brooks expression core */
/**
 * @file csf_pair_core.c
 * @brief Per-pair coupling core: @ref pair_eval_internal, the single-(bra,ket)
 *        matrix-element evaluator the whole @c csf/ public API funnels into.
 *
 * @c csf_pair.c calls @ref pair_eval_internal once per CSF pair (on a tiny
 * sub-DRT). This file is the DOWNSTREAM OWNER of the pair-coupling physics: it
 * runs the Brooks case catalog (@ref case01_pair .. @ref case14_pair), each a
 * tensor product walked by @c calc_TensorOp_average (in @c expr_tree.c), and
 * folds the resulting loop elements into the terms of ONE (bra,ket) coupling via
 * @ref add_expr_pair.
 *
 * It is the pair-evaluation twin of @c expr.c's full-expression path
 * (@c brooks_all / @c add_expr / @c expand_expr): same Brooks cases, same tensor
 * products and coefficients, but instead of emitting the whole expression it
 * filters loop elements down to the single (bra,ket) matrix element with two
 * pair-specific "D3" filters (see @ref add_expr_pair).
 *
 * Ownership of deferred questions, NOT restated here:
 *   - segment-factor values (the @c matrix_el / @c hmlkva evaluation reached
 *     through @c calc_TensorOp_average) are owned by @c expr_tree.c's todos;
 *   - the DRT arc-weight / canonical-index convention by @c drt.c
 *     (@c store_number_upper_arcs) and @c csf_lookup.c;
 *   - the @c n_oel one-/two-electron boundary and the @c ij pair-address meaning
 *     by @c expr.c (@c hexpr_expr_build / @c add_expr).
 * The <tt>\@todo domain:</tt> markers below cover only what this file itself
 * originates: the per-case Brooks coupling choices and the pair-extraction
 * filtering. The per-case questions overlap @c expr.c's @c caseNN todos (the two
 * paths share the Brooks catalog) and are cross-referenced as such.
 */
#include "csf_pair_core.h"
#include "drt.h"
#include "expr_internal.h"
#include "brooks_table.h"
#include "lk_table.h"
#include "dhash.h"
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

/* hexpr_init() (declared in hexpr.h via csf_pair_core.h) is idempotent;
 * pair_eval_internal calls it to ensure Wigner tables are ready before
 * calc_TensorOp_average. */

/* Operator codes (OP_UNITY..OP_ACT) are shared in expr_internal.h. */

/** @brief Seed the shared roots (@ref hexpr_roots_init) before any walk.
 *  The table driver reads coefficients from @c hexpr_root_val via
 *  @ref bc_fill_coef, so no file-static S2/S3 is needed. Called once at the
 *  top of @ref pair_eval_internal. */
static void init_constants(void) { hexpr_roots_init(); }

/** @brief Larger of two integers. @param a,b operands @return max(a,b). */
static int imax2(int a, int b) { return a > b ? a : b; }
/** @brief Smaller of two integers. @param a,b operands @return min(a,b). */
static int imin2(int a, int b) { return a < b ? a : b; }
/**
 * @brief Symmetric packed pair index @f$ \mathrm{max}(\mathrm{max}+1)/2+\mathrm{min} @f$.
 * @param a,b the two values (order-independent).
 * @return the packed index. Used to form the @c ij row address of a term; the
 *         meaning of @c ij as a CSF-pair address is owned by @c expr.c.
 */
static int triangle(int a, int b) { int mx=imax2(a,b),mn=imin2(a,b); return mx*(mx+1)/2+mn; }
/**
 * @brief Strictly-upper-triangular packed index @f$ \mathrm{max}(\mathrm{max}-1)/2+\mathrm{min} @f$.
 * @param p,q the two values (order-independent).
 * @return the packed index. Building block for the orbital @c seq_num (pqrs) in
 *         @ref add_expr_pair; the pqrs ordering convention is owned by @c expr.c.
 */
static int tri(int p, int q)      { int mx=imax2(p,q),mn=imin2(p,q); return mx*(mx-1)/2+mn; }

/**
 * @brief Build a ::tensor_product_t from positional literals (so the @c caseNN
 *        routines read as a compact table).
 * @param lenop number of active operator slots (1..4).
 * @param o0,o1,o2,o3 operator codes (@c OP_UNITY etc.) for slots 0..3.
 * @param t0,t1,t2,t3 the @c theta value for slots 0..3.
 * @return the initialized tensor product (delegates to @c tp_init). Mirrors the
 *         identical @c TP helper in @c expr.c.
 */
static tensor_product_t TP(int lenop,
                            int o0, int o1, int o2, int o3,
                            int t0, int t1, int t2, int t3)
{
    tensor_product_t tp;
    int op[4]    = {o0, o1, o2, o3};
    int theta[4] = {t0, t1, t2, t3};
    tp_init(&tp, op, theta, lenop);
    return tp;
}

#define CALC(drt, tp, origin) calc_TensorOp_average((drt), &(tp), (origin))

/* per-pair dedup entry */
typedef struct {
    int bottom, top, weight_g, weight_f, seq_num;
    double value;
} pair_term_t;

/** @brief Value-destructor passed to @c dhash_free: frees one heap
 *  ::pair_term_t. @param p the stored value pointer. */
static void free_pt(void *p) { free(p); }

/* Filter by D3, dedup by (ij_key, seq_num), expand for exact (bra, ket).
 * Each le_vecs[i] contributes with mul_coefs[i]; all share the same a[].
 *
 * Note: although the sub-DRT passed by csf_pair.c contains only nodes
 * reachable by the two input CSFs, the D3 filters below remain necessary.
 * Measured: the delta filter passes only ~3-7% of loop
 * elements; the path-weight filter passes only ~1.5-5%. The filters reject
 * the directionally-opposite half of every loop element and also reject
 * elements arising from cross-product paths (section6.2 of analysis.md) that the
 * 2-CSF sub-DRT admits through shared intermediate nodes. */
/**
 * @brief Fold the loop elements of one or more tensor products into the terms of
 *        a SINGLE (bra, ket) coupling, applying the two pair-extraction filters.
 * @param drt       the (sub-)DRT the loop elements were walked from.
 * @param le_vecs   array of @p nvec loop-element vectors (one per tensor product).
 * @param mul_coefs per-vector scalar multiplier; @c mul_coefs[i] scales
 *                  @c le_vecs[i].
 * @param nvec      number of vectors.
 * @param a         index map selecting which @c ijkl entries form the orbital key
 *                  (a pair when @p na<=2, a quartet when @p na>2).
 * @param na        number of orbital indices (2 or 4).
 * @param bra,ket   the canonical CSF indices whose coupling is being extracted.
 * @return an ::ot_vec_t of terms for this (bra,ket); caller owns @c .data.
 *
 * Role: the pair-evaluation analog of @c expr.c's @c add_expr + @c expand_expr,
 * specialized to ONE matrix element. It is the workhorse every @c caseNN_pair
 * calls. Data flow: accumulate over loop elements into a dedup hash keyed on
 * (ij_key, seq_num), summing coefficient-scaled values and dropping terms that
 * cancel below @c 1e-5; then expand each surviving term over the path-weight
 * pairs that actually reach @p bra, stamping the global @c triangle(bra,ket) as
 * @c ij.
 *
 * Two pair-specific "D3" filters (see the block comment above) reject loop
 * elements that do not connect exactly this bra and ket: (1) the @b delta filter
 * keeps only elements whose @c weight_g-weight_f equals @c bra-ket (the right
 * direction/offset); (2) the @b path-weight filter keeps only elements that have
 * some lower/upper path-weight pair summing to @c bra-weight_g (a walk actually
 * landing on @p bra). The expansion pass re-applies the same @c ==bra test so
 * only the matching paths emit terms.
 *
 * @note The @c seq_num (pqrs) construction and the @c n_oel split follow the
 *       conventions owned by @c expr.c's @c add_expr; the arc/path-weight reads
 *       follow @c drt.c / @c csf_lookup.c. Those are not re-questioned here.
 *
 * DOMAIN (per paper + Shavitt DRT): The (bra,ket) coupling depends only on the
 * bra and ket walks (per-orbital factorization, Sec. 2.4, cf. csf_pair.c), so the
 * loop elements contributing to <bra|H|ket> are exactly those whose two walks are
 * bra and ket. The 2-CSF sub-DRT also admits the cross-product pairs (bra-bra,
 * ket-ket) and the reversed (ket,bra); the D3 filters select the wanted pair:
 * (1) the delta test weight_g-weight_f == bra-ket keeps only the signed spin
 * offset of the target pair, rejecting the direction-reversed half; (2) the
 * landing test -- lower+upper path-weight summing to bra-weight_g -- requires the
 * walk to actually address bra in the Shavitt reverse-lexical system
 * (cf. csf_lookup.c), rejecting cross-product paths. Together they isolate the
 * (bra,ket) loop elements. See Sec. 2.4; addressing cf. csf_lookup.c.
 */
static ot_vec_t add_expr_pair(const hexpr_drt_t *drt,
                               le_vec_t **le_vecs, const double *mul_coefs, int nvec,
                               const int *a, int na,
                               int bra, int ket)
{
    int norb  = hexpr_drt_norb(drt);
    int n_oel = norb * (norb + 1) / 2;
    int delta = bra - ket;
    dhash_t *h = dhash_create();

    /* Outer: one input vector per tensor product, scaled by mul_coefs[vi].
     * Inner: every loop element produced for that product. */
    for (int vi = 0; vi < nvec; vi++) {
        le_vec_t *lev = le_vecs[vi];
        double mc = mul_coefs[vi];
        for (int j = 0; j < lev->size; j++) {
            loop_element_t *x = &lev->data[j];
            /* D3 filter 1 (delta): keep only elements whose g/f weight offset
             * matches the bra-ket direction; rejects the opposite half. */
            if (x->weight_g - x->weight_f != delta) continue;

            int nlw, nuw;
            const int *lw = hexpr_drt_lower_path_weights(drt, x->bottom, &nlw);
            const int *uw = hexpr_drt_upper_path_weights(drt, x->top,    &nuw);
            if (nlw == 0 || nuw == 0) continue;

            /* D3 filter 2 (path-weight): keep only elements that have some
             * lower+upper path-weight pair landing exactly on bra. */
            int target_sum = bra - x->weight_g;
            int passes = 0;
            for (int lx = 0; lx < nlw && !passes; lx++)
                for (int uy = 0; uy < nuw && !passes; uy++)
                    if (lw[lx] + uw[uy] == target_sum) passes = 1;
            if (!passes) continue;

            /* seq_num = packed orbital (pqrs) address (convention owned by
             * expr.c's add_expr): a pair index when na<=2, else a nested quartet
             * index shifted past the n_oel one-electron block. */
            int seq_num;
            if (na > 2) {
                int t0=x->ijkl[a[0]], t1=x->ijkl[a[1]];
                int t2=x->ijkl[a[2]], t3=x->ijkl[a[3]];
                seq_num = tri(tri(t0,t1), tri(t2,t3)) - 1 + n_oel;
            } else {
                seq_num = tri(x->ijkl[a[0]], x->ijkl[a[1]]) - 1;
            }

            /* Accumulate into the dedup hash keyed on (ij_key, seq_num); the
             * representative weight uses the [0] path-weight pair, matching
             * expr.c's add_expr. */
            double val = mc * x->value;
            int w = lw[0] + uw[0];
            int ij_key = triangle(w + x->weight_g, w + x->weight_f);
            dhash_key_t key = {(int64_t)ij_key, (int64_t)seq_num};

            pair_term_t *existing = (pair_term_t *)dhash_get(h, key);
            if (existing) {
                /* Same address: sum the contribution; evict if it cancels out. */
                existing->value += val;
                if (fabs(existing->value) <= 1e-5) {
                    free(existing);
                    dhash_delete(h, key);
                }
            } else if (fabs(val) > 1e-5) {
                pair_term_t *y = malloc(sizeof(pair_term_t));
                y->bottom   = x->bottom;
                y->top      = x->top;
                y->weight_g = x->weight_g;
                y->weight_f = x->weight_f;
                y->seq_num  = seq_num;
                y->value    = val;
                dhash_set(h, key, y);
            }
        }
    }

    ot_vec_t out = {NULL, 0, 0};
    int n = dhash_size(h);
    if (n > 0) {
        dhash_key_t *keys = malloc(n * sizeof(dhash_key_t));
        void **vals = malloc(n * sizeof(void *));
        dhash_to_array(h, keys, vals);
        int target_ij = triangle(bra, ket);
        /* Expansion pass: re-apply the path-weight==bra test to each surviving
         * term, emitting one output term (stamped with the global bra/ket ij)
         * per lower/upper path-weight pair that lands on bra. */
        for (int i = 0; i < n; i++) {
            pair_term_t *pt = (pair_term_t *)vals[i];
            int nlw, nuw;
            const int *lw = hexpr_drt_lower_path_weights(drt, pt->bottom, &nlw);
            const int *uw = hexpr_drt_upper_path_weights(drt, pt->top,    &nuw);
            for (int lx = 0; lx < nlw; lx++) {
                for (int uy = 0; uy < nuw; uy++) {
                    if (lw[lx] + pt->weight_g + uw[uy] == bra) {
                        out_term_t ot = {target_ij, pt->seq_num, pt->value};
                        ot_push(&out, &ot);
                    }
                }
            }
        }
        free(keys);
        free(vals);
    }
    dhash_free(h, free_pt);
    return out;
}

/**
 * @brief Move every term from @p src onto the end of @p all, then empty/free @p src.
 * @param all destination accumulator (keeps ownership of its buffer).
 * @param src source vector; consumed -- its @c data is freed and it is reset.
 * @note "Move", not copy: @p src is left empty so its buffer is never
 *       double-freed. Used by @ref pair_eval_internal and the dual-output cases.
 */
static void append_ot(ot_vec_t *all, ot_vec_t *src)
{
    for (int i = 0; i < src->size; i++) ot_push(all, &src->data[i]);
    free(src->data);
    src->data = NULL; src->size = src->cap = 0;
}

/**
 * @brief Evaluate the full coupling of one (bra, ket) CSF pair: run every Brooks
 *        case and concatenate their terms.
 * @param drt the (sub-)DRT addressing @p bra / @p ket (built by @c csf_pair.c).
 * @param bra,ket canonical CSF indices of the pair.
 * @return an ::ot_vec_t of all contributing terms; the CALLER owns @c .data.
 *
 * Role: the public entry point of this unit -- @c csf_pair.c's
 * @c pair_via_subdrt calls it once per pair. It ensures the Wigner tables are
 * built (@c hexpr_init, idempotent) and the shared roots are seeded
 * (@c hexpr_roots_init), then iterates @ref BROOKS_CASES, walks each enabled
 * row's ops, runs each term through @ref add_expr_pair, and @ref append_ot
 * "appends" the terms onto one growing result. NO final sort: @c ij is constant
 * (@c triangle(bra,ket)) for a pair, so emission order is the output order.
 *
 * @note The case set run here is the pair-path analog of @c expr.c's
 *       @c brooks_all but is NOT identical: this path includes the
 *       @c enable_pair rows -- so it RUNS case03 (which @c brooks_all omits via
 *       @c enable_whole==0) and, like @c brooks_all, skips the disabled case9.
 *
 * DOMAIN (per paper + code reference): The redundancy that reduces Brooks' 14
 * cases (Sec. 3.1, Table 1) to a generating set comes from the two-electron
 * integral's permutation symmetries and the ID(Psi) <= ID(Phi) convention
 * (Sec. 3.1). This per-pair path keeps case 03, unlike the full-expression
 * brooks_all: the full path visits ordered pairs under ID(Psi) <= ID(Phi), so one
 * of the two one-electron directions suffices and the other follows by transpose;
 * the per-pair path evaluates an arbitrary (bra,ket) directly and needs both, so
 * it runs case03 (cf. case03_pair). Case 09 is disabled in both paths. The paper
 * gives the redundancy principle (Sec. 3.1); the exact per-path case set is
 * recorded in the code's analysis (ruby-analysis.md / run_flag / spec D1), not in
 * the paper. See Sec. 3.1; case set per ruby-analysis.md.
 */
ot_vec_t pair_eval_internal(const hexpr_drt_t *drt, int bra, int ket)
{
    hexpr_init();
    init_constants();          /* seeds hexpr_root_val (shared roots) before any walk */
    hexpr_lk_ensure_init();    /* shared L(k) table must be initialized before the hot line */
    /* The pair path shares the SAME L(k) table as the whole path but stays in
     * OFF: it is NOT warmed, so the hot line fills-on-miss (uncounted). This is
     * correct -- the enable_pair rows (incl. case03) on per-pair sub-DRTs may need
     * keys the whole-path warm never filled. OFF fill-on-miss is numerically
     * identical to warm-fill (hmlkva is pure), so pair misses are expected, not a
     * failure. Only the WHOLE search pass asserts misses==0. The whole path always
     * ends in finish() (=> OFF), and the pair path never changes the phase, so the
     * table is always OFF on entry here; assert that invariant explicitly. */
    assert(hexpr_lk_phase() == LK_OFF);
    ot_vec_t all = {NULL, 0, 0};

    for (int c = 0; c < N_BROOKS_CASES; c++) {          /* cases */
        const brooks_case_t *bc = &BROOKS_CASES[c];
        if (!bc->enable_pair) continue;                /* skips only case9 (case03 runs) */

        /* Walk every op once into ev[]. */
        le_vec_t ev[BC_MAX_OPS];
        for (int o = 0; o < bc->nops; o++) {           /* ops */
            tensor_product_t tp = TP(bc->ops[o].lenop,
                bc->ops[o].chain[0], bc->ops[o].chain[1],
                bc->ops[o].chain[2], bc->ops[o].chain[3],
                bc->ops[o].theta[0], bc->ops[o].theta[1],
                bc->ops[o].theta[2], bc->ops[o].theta[3]);
            ev[o] = CALC(drt, tp, bc->origin);
        }

        /* Each term = one add_expr_pair over the ops it selects. */
        for (int t = 0; t < bc->nterms; t++) {         /* terms */
            const brooks_term_t *tm = &bc->terms[t];
            le_vec_t *sel[BC_MAX_OPS];
            for (int i = 0; i < tm->n_op_indices; i++)
                sel[i] = &ev[tm->op_indices[i]];
            double coef[BC_MAX_OPIDX];
            bc_fill_coef(tm, coef);
            ot_vec_t r = add_expr_pair(drt, sel, coef, tm->n_op_indices,
                                       tm->perm, tm->perm_len, bra, ket);
            append_ot(&all, &r);                       /* emission order preserved */
        }

        for (int o = 0; o < bc->nops; o++) free(ev[o].data);
    }
    return all;   /* NO final sort: ij is constant (triangle(bra,ket)) per pair,
                   * so emission order is the output order. */
}

