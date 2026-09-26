/* loopgen/expr.c -- BrooksCases, add_expr, expand_expr, hexpr_expr_t */
/**
 * @file expr.c
 * @brief Coupling-coefficient expression generator for the GUGA pipeline.
 *
 * This translation unit turns a fully-built Distinct Row Table (DRT) into a
 * flat, sorted list of energy-expression terms. The data flow is a three-stage
 * pipeline, each stage feeding the next:
 *
 *   1. The @c case01 .. @c case14 routines each describe one "Brooks case" as a
 *      tensor product of operator codes (see @ref OP_UNITY etc.) and hand it to
 *      @ref calc_TensorOp_average (implemented in expr_tree.c), which walks the
 *      DRT and returns the matching loop elements (::loop_element_t).
 *   2. @ref add_expr accumulates those loop elements into compressed terms
 *      (::comp_term_t), keyed by an (ij_key, seq_num) pair, summing coefficients
 *      and dropping anything below a numerical threshold.
 *   3. @ref expand_expr re-expands each compressed term over every lower/upper
 *      path-weight pair into the final output terms (::out_term_t).
 *
 * @ref brooks_all drives stages 1-3 for all enabled cases and sorts the result;
 * @ref hexpr_expr_build wraps that together with the CSF list into the opaque
 * ::hexpr_expr_t handle exposed by the public API.
 *
 * @note The integer "case" labels and operator orderings mirror the original
 *       Ruby reference (emol_const_ops.rb / BrooksCases); the in-body comments
 *       below describe only the data flow visible in this file. The GUGA /
 *       Paldus-Shavitt physics behind each coefficient is flagged with
 *       <tt>\@todo domain:</tt> markers for expert sign-off and is not asserted
 *       here. Those markers aggregate into Doxygen's Todo List page.
 */
#include "expr.h"
#include "expr_internal.h"
#include "brooks_table.h"
#include "fmt.h"
#include "dhash.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Operator codes (OP_UNITY..OP_ACT) are shared in expr_internal.h. */

/* ------------------------------------------------------------------ */
/* Math helpers                                                        */
/* ------------------------------------------------------------------ */
/** @brief Larger of two integers. @param a,b operands @return max(a,b). */
static int imax2(int a, int b) { return a > b ? a : b; }
/** @brief Smaller of two integers. @param a,b operands @return min(a,b). */
static int imin2(int a, int b) { return a < b ? a : b; }

/* triangle(a,b) = max(a,b)*(max(a,b)+1)/2 + min(a,b) */
/**
 * @brief Symmetric 0-based packed index of an unordered pair (a,b).
 * @param a,b the two integers (order does not matter).
 * @return @f$ \mathrm{max}(\mathrm{max}+1)/2 + \mathrm{min} @f$, where
 *         max/min are taken over {a,b}.
 * @note Used to fold the two path weights of a term into a single @c ij_key.
 */
static int triangle(int a, int b)
{
    int mx = imax2(a, b), mn = imin2(a, b);
    return mx * (mx + 1) / 2 + mn;
}

/* tri(p,q) = max(p,q)*(max(p,q)-1)/2 + min(p,q) */
/**
 * @brief Strictly-upper-triangular packed index of an unordered pair (p,q).
 * @param p,q the two integers (order does not matter).
 * @return @f$ \mathrm{max}(\mathrm{max}-1)/2 + \mathrm{min} @f$.
 * @note Differs from @ref triangle by the @c -1 (excludes the diagonal); used
 *       to build the orbital index @c seq_num inside @ref add_expr.
 */
static int tri(int p, int q)
{
    int mx = imax2(p, q), mn = imin2(p, q);
    return mx * (mx - 1) / 2 + mn;
}

/* ------------------------------------------------------------------ */
/* add_expr                                                            */
/* ------------------------------------------------------------------ */
/* Returns a vector of comp_term_t entries (sorted by [ij_key, seq_num]).
 * Modifies seq_num on the loop elements in le_vecs (like Ruby mutates x). */

/**
 * @brief qsort/stable-sort comparator ordering ::comp_term_t by (ij_key, seq_num).
 * @param a,b pointers to ::comp_term_t.
 * @return negative/zero/positive following the (ij_key, then seq_num) ordering.
 */
static int ct_cmp_key(const void *a, const void *b)
{
    const comp_term_t *ca = (const comp_term_t *)a;
    const comp_term_t *cb = (const comp_term_t *)b;
    if (ca->ij_key != cb->ij_key) return ca->ij_key - cb->ij_key;
    return ca->seq_num - cb->seq_num;
}

/* Used only for the final stable sort across all cases: sort by ij_key only,
 * preserving per-case (ij_key, seq_num) order within equal ij_key groups.
 * Matches Ruby: a.sort{ |x,y| x[0][0] - y[0][0] } */
/**
 * @brief qsort/stable-sort comparator ordering ::comp_term_t by @c ij_key only.
 * @param a,b pointers to ::comp_term_t.
 * @return ca->ij_key - cb->ij_key.
 * @note Relies on a *stable* sort to preserve per-case (ij_key, seq_num) order
 *       within an equal-ij_key group; see the call site in @ref brooks_all.
 */
static int ct_cmp_ij_only(const void *a, const void *b)
{
    return ((const comp_term_t *)a)->ij_key - ((const comp_term_t *)b)->ij_key;
}

/* Allocated comp_term_t used as dhash value */
/** @brief Value-destructor passed to @c dhash_free: frees one heap ::comp_term_t.
 *  @param p the stored value pointer (a @c comp_term_t* upcast to @c void*). */
static void free_ct(void *p) { free(p); }

/**
 * @brief Accumulate the loop elements of one or more tensor products into a
 *        sorted, de-duplicated vector of compressed expression terms.
 *
 * Stage 2 of the pipeline (see @ref expr.c "file overview"). Each input
 * ::loop_element_t carries a DRT walk endpoint (@c bottom / @c top), a pair of
 * path-weight offsets (@c weight_g / @c weight_f), orbital labels @c ijkl, and a
 * numeric @c value. This routine collapses elements that land on the same
 * (ij_key, seq_num) address by summing their (coefficient-scaled) values, and
 * discards any accumulated term whose magnitude falls to/below @c 1e-5.
 *
 * @param drt     fully-built DRT the loop elements were walked from.
 * @param le_vecs array of @p nvec loop-element vectors (one per tensor product).
 * @param nvec    number of vectors in @p le_vecs.
 * @param coef    per-vector scalar multiplier; @c coef[i] scales every element
 *                of @c le_vecs[i].
 * @param a       index map selecting which of a loop element's @c ijkl[] entries
 *                form the orbital key: a pair when @p na <= 2, a quartet when
 *                @p na > 2.
 * @param na      number of orbital indices the key is built from (2 or 4).
 * @return a freshly-allocated ::ct_vec_t sorted by (ij_key, seq_num); the caller
 *         owns @c out.data. Empty (NULL data) when nothing survived the threshold.
 *
 * @note Side effect: writes each visited element's @c seq_num field, mirroring
 *       the in-place mutation the Ruby reference performs on its loop elements.
 */
static ct_vec_t add_expr(const hexpr_drt_t *drt,
                          le_vec_t *le_vecs, int nvec,
                          const double *coef, const int *a, int na)
{
    int norb  = hexpr_drt_norb(drt);
    int n_oel = norb * (norb + 1) / 2;

    /* Hash table keyed by (ij_key, seq_num); values are heap comp_term_t.
     * Using a hash here is what makes the accumulation O(elements) rather than
     * an O(n^2) search for an existing term at the same address. */
    dhash_t *h = dhash_create();

    /* Outer loop: one input vector per tensor product, sharing scalar coef[i].
     * Inner loop: every loop element produced by the DRT walk for that product. */
    for (int i = 0; i < nvec; i++) {
        le_vec_t *lev = &le_vecs[i];
        for (int j = 0; j < lev->size; j++) {
            loop_element_t *x = &lev->data[j];

            /* Resolve the element's DRT endpoints into path-weight lists. An
             * endpoint with no reachable paths contributes nothing -- skip it
             * before it can corrupt the index arithmetic below. */
            int nlw, nuw;
            const int *lw = hexpr_drt_lower_path_weights(drt, x->bottom, &nlw);
            const int *uw = hexpr_drt_upper_path_weights(drt, x->top,    &nuw);
            if (nlw == 0 || nuw == 0) continue;

            /* Representative weight: first lower + first upper path weight. The
             * full fan-out over all path-weight pairs is deferred to
             * @ref expand_expr; here only the [0] pair seeds the dedupe key. */
            int w = lw[0] + uw[0];

            /* compute seq_num from ijkl[a[...]] -- mutates x->seq_num */
            /* seq_num is the packed orbital (pqrs) address of this term: a
             * one-electron pair index when na<=2, otherwise a two-electron
             * quartet index nested as tri(tri(.,.),tri(.,.)) and shifted past
             * the n_oel one-electron block. */
            int seq_num;
            if (na > 2) {
                int t0 = x->ijkl[a[0]], t1 = x->ijkl[a[1]];
                int t2 = x->ijkl[a[2]], t3 = x->ijkl[a[3]];
                seq_num = tri(tri(t0,t1), tri(t2,t3)) - 1 + n_oel;
            } else {
                int t0 = x->ijkl[a[0]], t1 = x->ijkl[a[1]];
                seq_num = tri(t0, t1) - 1;
            }
            x->seq_num = seq_num;  /* mutate like Ruby */

            /* ij_key folds the two shifted path weights (w + bra/ket offsets)
             * into one symmetric packed index -- the "row" address of the term.
             * Together (ij_key, seq_num) is the full coordinate used for dedupe. */
            int ij_key = triangle(w + x->weight_g, w + x->weight_f);
            dhash_key_t key = {(int64_t)ij_key, (int64_t)seq_num};

            /** DOMAIN (per paper + Shavitt DRT): ij_key is the (bra,ket) CSF-pair
             *  address -- triangle() packs the two shifted path weights (each a
             *  Shavitt reverse-lexical CSF index, cf. csf_lookup.c) into one
             *  symmetric row index. seq_num is the integral label: the
             *  one-electron block uses tri(i,j), the two-electron block
             *  tri(tri..)+n_oel, i.e. the (pq|rs) index into the MO integral
             *  array. Together they form the output tuple
             *  (ID(Psi),ID(Phi),int_type,orbital indices) of the emission format
             *  (Sec. 4.3); the paper specifies this layout, the pair-address
             *  packing follows Shavitt. */
            comp_term_t *existing = (comp_term_t *)dhash_get(h, key);
            if (existing) {
                /* Term already at this address: fold in the new contribution.
                 * If the running sum cancels below threshold, evict the entry
                 * entirely so it never reaches the output. */
                existing->value += coef[i] * x->value;
                /* apply threshold */
                if (fabs(existing->value) > 1e-5)
                    dhash_set(h, key, existing);
                else {
                    free(existing);
                    dhash_delete(h, key);
                }
            } else {
                /* First contribution at this address: allocate, copy the
                 * element's endpoints/weights, and insert only if it already
                 * clears the threshold (a single tiny term is never stored). */
                comp_term_t *y = malloc(sizeof(comp_term_t));
                y->ij_key   = ij_key;
                y->seq_num  = seq_num;
                y->bottom   = x->bottom;
                y->top      = x->top;
                y->weight_g = x->weight_g;
                y->weight_f = x->weight_f;
                y->value    = coef[i] * x->value;
                /* threshold: |y->value| must be > 1e-5 */
                if (fabs(y->value) > 1e-5)
                    dhash_set(h, key, y);
                else
                    free(y);
            }
        }
    }

    /* Drain phase: the hash now holds exactly the surviving terms. Snapshot
     * its live entries into a flat array, sort into canonical (ij_key, seq_num)
     * order, then release the hash (values freed via @ref free_ct). */
    /* collect, sort by (ij_key, seq_num), and free hash */
    int n = dhash_size(h);
    ct_vec_t out = {NULL, 0, 0};
    if (n > 0) {
        dhash_key_t *keys = malloc(n * sizeof(dhash_key_t));
        void **vals       = malloc(n * sizeof(void *));
        dhash_to_array(h, keys, vals);

        out.data = malloc(n * sizeof(comp_term_t));
        out.size = n;
        out.cap  = n;
        for (int i = 0; i < n; i++)
            out.data[i] = *(comp_term_t *)vals[i];

        /* sort by [ij_key, seq_num] -- stable; ties impossible since keys unique */
        hexpr_stable_sort(out.data, out.size, sizeof(comp_term_t), ct_cmp_key);

        free(keys);
        free(vals);
    }
    dhash_free(h, free_ct);
    return out;
}

/* ------------------------------------------------------------------ */
/* expand_expr                                                         */
/* ------------------------------------------------------------------ */
/* Input: sorted ct_vec_t.  Output: sorted ot_vec_t.                  */

/**
 * @brief qsort/stable-sort comparator ordering ::out_term_t by @c ij.
 * @param a,b pointers to ::out_term_t.
 * @return a->ij - b->ij.
 */
static int ot_cmp_ij(const void *a, const void *b)
{
    return ((const out_term_t *)a)->ij - ((const out_term_t *)b)->ij;
}

/**
 * @brief Expand each compressed term over every lower/upper path-weight pair
 *        into the final flat list of output terms.
 *
 * Stage 3 of the pipeline. @ref add_expr deliberately keyed terms on only the
 * first ([0]) path-weight pair; this routine restores the full set: for one
 * ::comp_term_t it emits one ::out_term_t per (lower_weight, upper_weight)
 * combination, recomputing the row index @c ij for each pair while carrying the
 * shared @c pqrs (= @c seq_num) and @c coef across all of them.
 *
 * @param drt fully-built DRT (source of the path-weight lists).
 * @param v   sorted compressed terms from @ref add_expr / @ref brooks_all.
 * @return a freshly-allocated ::ot_vec_t sorted by @c ij; caller owns @c out.data.
 *
 * @note The nested @c lx / @c uy loops are a Cartesian product over the two
 *       weight lists -- term count grows as @c nlw*nuw per input term. Endpoints
 *       with no reachable paths are skipped (same guard as @ref add_expr).
 *
 * DOMAIN (Shavitt DRT; addressing used by the paper, not re-derived there):
 * Each (lower_weight, upper_weight) pair maps to a separate ij row because the
 * two path weights address distinct CSFs: the lower weight indexes a walk below
 * the term's bottom node, the upper weight a walk above its top node, and their
 * combination is one (bra,ket) CSF pair in the Shavitt reverse-lexical system
 * (cf. csf_lookup.c). The Cartesian product over the two weight lists therefore
 * enumerates every CSF pair the loop contributes to, one ij row each. This is
 * the Shavitt path-weight addressing; the paper builds the DRT (Sec. 4.1)
 * without restating the correspondence.
 */
static ot_vec_t expand_expr(const hexpr_drt_t *drt, const ct_vec_t *v)
{
    ot_vec_t out = {NULL, 0, 0};

    for (int i = 0; i < v->size; i++) {
        const comp_term_t *t = &v->data[i];
        int nlw, nuw;
        const int *lw = hexpr_drt_lower_path_weights(drt, t->bottom, &nlw);
        const int *uw = hexpr_drt_upper_path_weights(drt, t->top,    &nuw);
        if (nlw == 0 || nuw == 0) continue;
        /* Cartesian product: every lower path weight against every upper one.
         * Each pair yields one concrete row index ij; pqrs and coef are
         * invariant across the expansion of a single compressed term. */
        for (int lx = 0; lx < nlw; lx++) {
            for (int uy = 0; uy < nuw; uy++) {
                out_term_t ot;
                ot.ij   = triangle(lw[lx] + uw[uy] + t->weight_g,
                                   lw[lx] + uw[uy] + t->weight_f);
                ot.pqrs = t->seq_num;
                ot.coef = t->value;
                ot_push(&out, &ot);
            }
        }
    }

    hexpr_stable_sort(out.data, out.size, sizeof(out_term_t), ot_cmp_ij);
    return out;
}

/* ------------------------------------------------------------------ */
/* BrooksCases helpers                                                 */
/* ------------------------------------------------------------------ */

/* Append sorted ct_vec result to a growing flat ct_vec */
/**
 * @brief Move every term from one ::ct_vec_t onto the end of an accumulator,
 *        growing the accumulator by doubling, then empty/free the source.
 * @param all destination accumulator (ownership of its buffer stays with it).
 * @param src source vector; consumed -- its @c data is freed and it is reset.
 * @note "Move", not copy: @p src is left empty so its buffer is never
 *       double-freed by the caller.
 */
static void append_ct(ct_vec_t *all, ct_vec_t *src)
{
    for (int i = 0; i < src->size; i++) {
        if (all->size == all->cap) {
            all->cap = all->cap ? all->cap * 2 : 256;
            all->data = realloc(all->data, all->cap * sizeof(comp_term_t));
        }
        all->data[all->size++] = src->data[i];
    }
    free(src->data);
    src->data = NULL; src->size = src->cap = 0;
}

/* Allocate a tensor product from literal initializers */
/**
 * @brief Build a ::tensor_product_t from positional literals (a by-value helper
 *        so the @c caseNN routines read as a compact table).
 * @param lenop number of active operator slots (1..4).
 * @param o0,o1,o2,o3 operator codes (see @ref OP_UNITY etc.) for slots 0..3.
 * @param t0,t1,t2,t3 the @c theta value for slots 0..3.
 * @return the initialized tensor product (delegates to @c tp_init, which also
 *         derives the @c theta_prev chain).
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

/* ------------------------------------------------------------------ */
/* BrooksCases: case01..case14                                         */
/* ------------------------------------------------------------------ */

/** @brief Seed the shared roots (@ref hexpr_roots_init) before any walk.
 *  The table driver reads coefficients from @c hexpr_root_val via
 *  @ref bc_fill_coef, so no file-static S2/S3 is needed. Called once at the
 *  top of @ref brooks_all. */
static void bc_init_constants(void)
{
    hexpr_roots_init();
}

/* Macro to run one TensorProduct and collect loop elements */
#define CALC(drt, tp, origin) calc_TensorOp_average((drt), &(tp), (origin))


/* ------------------------------------------------------------------ */
/* warm_lk_table -- generate pass (fill the L(k) table)                */
/* ------------------------------------------------------------------ */
/**
 * @brief Fill the shared L(k) table by walking exactly the ops the search pass
 *        (@ref brooks_all's case loop) walks -- ops only, no add_expr.
 *
 * Port of RubyRef/emol_expression.rb warm_lk_table (:853-860). With the table in
 * GENERATE phase the hot line (@c lk_value in expr_tree.c) computes each
 * first-seen key's L(k) via hmlkva and inserts it. Because GENERATE computes the
 * real L(k) on a miss, the running product / threshold pruning is identical to
 * the search pass, so the tree shape -- and hence the key set -- matches exactly:
 * the search pass then finds every key present (misses==0).
 */
static void warm_lk_table(const hexpr_drt_t *drt)
{
    for (int c = 0; c < N_BROOKS_CASES; c++) {
        const brooks_case_t *bc = &BROOKS_CASES[c];
        if (!bc->enable_whole) continue;               /* same rows the case loop walks */
        for (int o = 0; o < bc->nops; o++) {           /* ops only -- no terms/add_expr */
            tensor_product_t tp = TP(bc->ops[o].lenop,
                bc->ops[o].chain[0], bc->ops[o].chain[1],
                bc->ops[o].chain[2], bc->ops[o].chain[3],
                bc->ops[o].theta[0], bc->ops[o].theta[1],
                bc->ops[o].theta[2], bc->ops[o].theta[3]);
            le_vec_t ev = CALC(drt, tp, bc->origin);   /* hot line fills LK (GENERATE) */
            free(ev.data);                             /* discard the walk result */
        }
    }
}

/* ------------------------------------------------------------------ */
/* BrooksCases::all                                                    */
/* ------------------------------------------------------------------ */
/* Returns expanded expression sorted by ij (stable). */

/**
 * @brief Drive every enabled Brooks case and assemble the final expression.
 *
 * Table-driven whole-matrix path: iterate @ref BROOKS_CASES, skip rows with
 * @c enable_whole==0 (case03 and the disabled case9), walk each of a row's ops
 * once with @c CALC, run each term through @ref add_expr, and append the result
 * onto one growing @c flat buffer. Then (unchanged) stable-sort @c flat by
 * @c ij_key alone and hand it to @ref expand_expr for the path-weight fan-out.
 *
 * @param drt fully-built DRT to evaluate.
 * @return the expanded, ij-sorted ::ot_vec_t; caller owns @c .data.
 *
 * The enabled-row order reproduces the previous hand-written emission order
 * exactly: case01, case02, case1..case8, case10..case14. Each term maps 1:1 to a
 * former @c add_expr call; coefficients are materialized from the row's
 * @c {mult,root} pairs via @ref bc_fill_coef through the shared @c hexpr_root_val.
 *
 * @note The final sort MUST be stable: @ref ct_cmp_ij_only compares @c ij_key
 *       only, relying on stability to keep each case's intra-group order.
 */
static ot_vec_t brooks_all(const hexpr_drt_t *drt)
{
    bc_init_constants();   /* seeds hexpr_root_val (shared roots) before any walk */

    /* Two-phase L(k) (Kamuy): pass 1 fills the table, pass 2 reads it.
     *   generate + warm : walk the enable_whole ops, filling LK via hmlkva.
     *   search          : the case loop below walks again; the hot line is a pure
     *                     table lookup -- we assert misses/wig/mel == 0 afterward.
     * WHOLE is STRICT (warm covers exactly what search reads => misses must be 0).
     * The PAIR path (csf_pair_core.c) is deliberately NOT warmed: it runs the
     * enable_pair rows (incl. case03) on per-pair sub-DRTs and legitimately needs
     * keys warm never filled, so it stays in OFF (fill-on-miss, uncounted) -- that
     * is correct, not a failure. Do not force the pair path through warm. */
    hexpr_lk_begin_generation();
    warm_lk_table(drt);
    hexpr_lk_begin_search();

    ct_vec_t flat = {NULL, 0, 0};

    for (int c = 0; c < N_BROOKS_CASES; c++) {          /* cases */
        const brooks_case_t *bc = &BROOKS_CASES[c];
        if (!bc->enable_whole) continue;               /* skips case03 and case9 */

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

        /* Each term = one add_expr over the ops it selects. */
        for (int t = 0; t < bc->nterms; t++) {         /* terms */
            const brooks_term_t *tm = &bc->terms[t];
            le_vec_t sel[BC_MAX_OPS];
            for (int i = 0; i < tm->n_op_indices; i++)
                sel[i] = ev[tm->op_indices[i]];
            double coef[BC_MAX_OPIDX];
            bc_fill_coef(tm, coef);
            ct_vec_t r = add_expr(drt, sel, tm->n_op_indices,
                                  coef, tm->perm, tm->perm_len);
            append_ct(&flat, &r);                      /* emission order preserved */
        }

        for (int o = 0; o < bc->nops; o++) free(ev[o].data);
    }

    hexpr_lk_finish();   /* search pass done; leave the shared table in OFF */

    /* final stable sort by ij_key only -- matches Ruby: a.sort{ |x,y| x[0][0] - y[0][0] }
     * within equal ij_key groups the per-case (ij_key,seq_num) order is preserved */
    hexpr_stable_sort(flat.data, flat.size, sizeof(comp_term_t), ct_cmp_ij_only);

    ot_vec_t out = expand_expr(drt, &flat);
    free(flat.data);
    return out;
}

/* ------------------------------------------------------------------ */
/* hexpr_expr_t                                                        */
/* ------------------------------------------------------------------ */

struct hexpr_expr {
    int     norb;
    int     ncsf;
    char   *csf_data;   /* norb * ncsf chars (NUL-terminated rows) */
    ot_vec_t terms;     /* full expression, sorted by ij */
    int     n_oel;      /* norb*(norb+1)/2: 1-electron threshold */
    int     n_one_terms; /* count of terms with pqrs < n_oel */
};

/**
 * @brief Build the opaque expression handle for a DRT: CSF list + sorted terms.
 *
 * Public entry point that packages the whole pipeline. It snapshots @c norb and
 * the one-electron block size @c n_oel, materializes the CSF strings via
 * @ref hexpr_drt_mk_csf_list into a single packed @c csf_data buffer, runs
 * @ref brooks_all to produce the sorted output terms, and precomputes how many
 * of those terms are one-electron (@c pqrs < @c n_oel) so the "one" views can be
 * served without rescanning.
 *
 * @param drt fully-built DRT to evaluate.
 * @return a heap ::hexpr_expr_t; release with @ref hexpr_expr_destroy (the
 *         name declared in @c hexpr.h); @ref hexpr_expr_free is the internal
 *         implementation it forwards to.
 *
 * @note Each CSF row is stored with stride @c norb+1 (the trailing NUL is kept),
 *       so row @c i lives at @c csf_data + i*(norb+1).
 */
hexpr_expr_t *hexpr_expr_build(const hexpr_drt_t *drt)
{
    hexpr_expr_t *e = calloc(1, sizeof(hexpr_expr_t));
    e->norb  = hexpr_drt_norb(drt);
    e->n_oel = e->norb * (e->norb + 1) / 2;

    /* Materialize the CSF list, then pack each NUL-terminated row contiguously
     * into csf_data (stride norb+1) and free the per-row temporaries. */
    /* CSF list */
    char **csfs = NULL;
    e->ncsf = hexpr_drt_mk_csf_list(drt, &csfs);
    e->csf_data = malloc((size_t)(e->ncsf) * (e->norb + 1));
    for (int i = 0; i < e->ncsf; i++) {
        memcpy(e->csf_data + i * (e->norb + 1), csfs[i], e->norb + 1);
        free(csfs[i]);
    }
    free(csfs);

    /* Run the full Brooks pipeline to obtain the ij-sorted output terms. */
    /* build expression */
    e->terms = brooks_all(drt);

    /* One pass to cache the one-electron term count (pqrs below the n_oel
     * boundary), so hexpr_expr_nterms_one / show_one need no rescan. */
    /* precompute one-electron term count */
    int n_one = 0;
    for (int i = 0; i < e->terms.size; i++)
        if (e->terms.data[i].pqrs < e->n_oel) n_one++;
    e->n_one_terms = n_one;

    return e;
}

/**
 * @brief Free an expression handle and everything it owns (CSF buffer, terms).
 * @param e handle to release; @c NULL is a no-op.
 */
void hexpr_expr_free(hexpr_expr_t *e)
{
    if (!e) return;
    free(e->csf_data);
    free(e->terms.data);
    free(e);
}

/** @brief Public alias for @ref hexpr_expr_free -- the name declared in
 *  @c hexpr.h and the only one the bindings call.
 *  @param e handle to release. */
void hexpr_expr_destroy(hexpr_expr_t *e) { hexpr_expr_free(e); }

/** @brief Total number of expression terms (one- and two-electron).
 *  @param e expression handle. @return term count as @c int64_t. */
int64_t hexpr_expr_nterms_full(const hexpr_expr_t *e)
{
    return (int64_t)e->terms.size;
}

/** @brief Number of one-electron terms (precomputed in @ref hexpr_expr_build).
 *  @param e expression handle. @return one-electron term count as @c int64_t. */
int64_t hexpr_expr_nterms_one(const hexpr_expr_t *e)
{
    return (int64_t)e->n_one_terms;
}

/**
 * @brief Copy expression terms into caller-provided parallel output arrays.
 * @param e expression handle.
 * @param which ::HEXPR_EXPR_FULL = all terms; ::HEXPR_EXPR_ONE = one-electron
 *        terms only (@c pqrs < n_oel).
 * @param[out] out_ij   per-term row index (0-based @c ij).
 * @param[out] out_pqrs per-term orbital index (@c pqrs / @c seq_num).
 * @param[out] out_coef per-term coefficient.
 * @return ::HEXPR_OK, or ::HEXPR_ERR_INVALID_ARG if any pointer is NULL.
 * @note The caller must size each array to at least @ref hexpr_expr_nterms_full
 *       (HEXPR_EXPR_FULL) or @ref hexpr_expr_nterms_one (HEXPR_EXPR_ONE)
 *       entries.
 */
hexpr_status_t hexpr_expr_read(const hexpr_expr_t *e,
                               hexpr_expr_which_t which,
                               int32_t *out_ij, int32_t *out_pqrs,
                               double *out_coef)
{
    if (!e || !out_ij || !out_pqrs || !out_coef) return HEXPR_ERR_INVALID_ARG;
    int k = 0;
    for (int i = 0; i < e->terms.size; i++) {
        const out_term_t *t = &e->terms.data[i];
        if (which == HEXPR_EXPR_ONE && t->pqrs >= e->n_oel) continue;
        out_ij[k]   = (int32_t)t->ij;
        out_pqrs[k] = (int32_t)t->pqrs;
        out_coef[k] = t->coef;
        k++;
    }
    return HEXPR_OK;
}

/* ------------------------------------------------------------------ */
/* Show helpers                                                        */
/* ------------------------------------------------------------------ */

/**
 * @brief Print the shared report preamble: the CSF list and the term-table
 *        column header.
 * @param e expression handle (source of the CSF rows).
 * @param f destination stream.
 * @return 0 on success, -1 on any @c fprintf I/O failure.
 */
static int write_header(const hexpr_expr_t *e, FILE *f)
{
    if (fprintf(f, " CSF List \n") < 0) return -1;
    for (int i = 0; i < e->ncsf; i++) {
        hexpr_print_csf(f, e->csf_data + i * (e->norb + 1), e->norb);
    }
    if (fprintf(f, " Energy Expression \n") < 0) return -1;
    if (fprintf(f, "      IJ       pqrs     Coef \n") < 0) return -1;
    return 0;
}

/**
 * @brief Print the full expression (header + every term) in fixed-column form.
 * @param e expression handle.
 * @param f destination stream.
 * @return ::HEXPR_OK, or ::HEXPR_ERR_IO on any write failure.
 * @note DELIBERATE dump convention: @c ij is printed 1-based (@c t->ij + 1) to
 *       match the Ruby reference corpus reference/<case>.expr.txt, while @c pqrs is
 *       printed 0-based (raw). The binary API (@ref hexpr_expr_read) emits BOTH
 *       0-based -- this dump-ij-1 / pqrs-0 asymmetry is intentional, not a bug.
 */
hexpr_status_t hexpr_expr_show_full(const hexpr_expr_t *e, FILE *f)
{
    if (write_header(e, f) < 0) return HEXPR_ERR_IO;
    for (int i = 0; i < e->terms.size; i++) {
        const out_term_t *t = &e->terms.data[i];
        /* t->ij + 1: deliberate 1-based dump ij (see @note); pqrs stays 0-based. */
        if (fprintf(f, " %7d    %10d    %24.16e \n",
                    t->ij + 1, t->pqrs, t->coef) < 0)
            return HEXPR_ERR_IO;
    }
    return HEXPR_OK;
}

/**
 * @brief Print only the one-electron terms (@c pqrs < n_oel), header included.
 * @param e expression handle.
 * @param f destination stream.
 * @return ::HEXPR_OK, or ::HEXPR_ERR_IO on any write failure.
 * @note Same format as @ref hexpr_expr_show_full, filtered to one-electron terms;
 *       same deliberate dump-ij-1 / pqrs-0 convention.
 */
hexpr_status_t hexpr_expr_show_one(const hexpr_expr_t *e, FILE *f)
{
    if (write_header(e, f) < 0) return HEXPR_ERR_IO;
    for (int i = 0; i < e->terms.size; i++) {
        const out_term_t *t = &e->terms.data[i];
        if (t->pqrs < e->n_oel) {
            /* t->ij + 1: deliberate 1-based dump ij (see @note); pqrs stays 0-based. */
            if (fprintf(f, " %7d    %10d    %24.16e \n",
                        t->ij + 1, t->pqrs, t->coef) < 0)
                return HEXPR_ERR_IO;
        }
    }
    return HEXPR_OK;
}
