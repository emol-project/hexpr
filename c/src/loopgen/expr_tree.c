/* loopgen/expr_tree.c -- hmlkva, matrix_el, tree_search, calc_TensorOp_average */
/**
 * @file expr_tree.c
 * @brief Loop-element generator: walks a DRT under a tensor product to produce
 *        the ::loop_element_t list that @c expr.c reduces into coupling terms.
 *
 * This unit closes the @c loopgen/ cross-link triangle: @c expr.c's @c caseNN
 * routines call @ref calc_TensorOp_average with a ::tensor_product_t, this code
 * walks the @ref hexpr_drt_t graph built in @c drt.c (reading it through the
 * @c hexpr_drt_* accessors), and the resulting loop elements flow back to
 * @c expr.c's @c add_expr.
 *
 * The data flow, bottom up:
 *   - @ref matrix_el and @ref hmlkva evaluate a single segment's numeric factor
 *     for a given walk-step pair and operator (the latter folding in a Wigner
 *     9j from @c wigner.h).
 *   - @ref tree_search recurses down the DRT, threading the per-level walk state
 *     (node indices, path weights, the running matrix-element product, the
 *     orbital-index stack) and emitting a ::loop_element_t each time a loop
 *     closes; a small record hash dedupes equal loops.
 *   - @ref calc_TensorOp_average drives that recursion from every admissible
 *     start node and returns the sorted, de-duplicated loop-element vector.
 *
 * @note In-body comments describe only the data flow visible in the code. The
 *       GUGA / Paldus-Shavitt-Racah intent -- the segment-factor values, the 9j
 *       coupling, why a closed walk contributes its product -- is deferred to
 *       <tt>\@todo domain:</tt> markers rather than asserted here.
 */
#include "expr_internal.h"
#include "lk_table.h"
#include "wigner.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Constants from emol_const_ops.rb                                    */
/* ------------------------------------------------------------------ */
/* Walk step codes (EMPTY..FULL) and operator codes (OP_UNITY..OP_ACT) are
 * now shared in expr_internal.h. */

/* The single sqrt(2)/sqrt(3) path for the whole expression generator.
 * hexpr_root_val[R1]=1.0 is exact; [R2]/[R3] are seeded here. Everything that
 * needs a root (matrix_el below; the table drivers' coef materialization) reads
 * this array, so there is exactly one place a root value can come from. */
double hexpr_root_val[4] = { 0.0, 1.0, 0.0, 0.0 };
void hexpr_roots_init(void)
{
    hexpr_root_val[R2] = sqrt(2.0);
    hexpr_root_val[R3] = sqrt(3.0);
}

/* ------------------------------------------------------------------ */
/* two-phase L(k) table (shared by the whole and pair paths)          */
/* ------------------------------------------------------------------ */
/* ONE process-lifetime table. The whole path fills it (generate/warm) then
 * reads it (search); the pair path leaves it in OFF and fills-on-miss. Because
 * L(k)'s key (gf,ff,hf,psi,phi,theta,theta_prev) is DRT-independent, the same
 * warmed value serves every DRT / sub-DRT, so one shared table is correct.
 * matrix_el and hmlkva read LK.phase to bump the SEARCH-phase call counters. */
static lk_table_t LK;
static int LK_inited = 0;

void hexpr_lk_ensure_init(void)   { if (!LK_inited) { lk_clear(&LK); LK_inited = 1; } }
void hexpr_lk_begin_generation(void) { hexpr_lk_ensure_init(); LK.phase = LK_GENERATE; }
void hexpr_lk_begin_search(void)  { LK.phase = LK_SEARCH; LK.miss = 0;
                                    LK.search_wig_calls = 0; LK.search_mel_calls = 0; }
void hexpr_lk_finish(void)        { LK.phase = LK_OFF; }
int  hexpr_lk_phase(void)         { return LK.phase; }
long hexpr_lk_search_misses(void)    { return LK.miss; }
long hexpr_lk_search_wig_calls(void) { return LK.search_wig_calls; }
long hexpr_lk_search_mel_calls(void) { return LK.search_mel_calls; }
long hexpr_lk_off_fills(void)        { return LK.off_fills; }
long hexpr_lk_table_size(void)       { return LK.size; }

/* REPR_STATE[gf]: (2j+1) rep; indices 0..4 */
static const int REPR_STATE[5] = {0, 1, 2, 2, 1};

/* REPR_OP[hf]: (2j+1) rank of operator; indices 0..13 */
static const int REPR_OP[14] = {0, 1, 2, 2, 1, 1, 3, 1, 2, 2, 1, 1, 3, 1};

/* NUM_CA[op]: creation+annihilation count; indices 0..13 */
static const int NUM_CA[14] = {0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 4, 1, 2, 2};

#define THRESH_ZERO 1.0e-7

/* Seg_walk[op][i] = {gf, ff} pairs; op index 0..10 */
static const int UNITY_W[][2] = {{1,1},{2,2},{3,3},{4,4},{2,3},{3,2}};
static const int A_W[][2]     = {{1,2},{1,3},{2,4},{3,4}};
static const int C_W[][2]     = {{2,1},{3,1},{4,2},{4,3}};
static const int AAS_W[][2]   = {{1,4}};
static const int CAS_W[][2]   = {{2,2},{3,3},{4,4},{2,3},{3,2}};
static const int CAT_W[][2]   = {{2,2},{3,3},{2,3},{3,2}};
static const int CCS_W[][2]   = {{4,1}};
static const int CAA_W[][2]   = {{2,4},{3,4}};
static const int CCA_W[][2]   = {{4,2},{4,3}};
static const int CCAA_W[][2]  = {{4,4}};

typedef struct { const int (*w)[2]; int n; } wset_t;
static const wset_t SEG_WALK[11] = {
    {NULL, 0},
    {UNITY_W, 6},
    {A_W,     4},
    {C_W,     4},
    {AAS_W,   1},
    {CAS_W,   5},
    {CAT_W,   4},
    {CCS_W,   1},
    {CAA_W,   2},
    {CCA_W,   2},
    {CCAA_W,  1},
};

/* ------------------------------------------------------------------ */
/* matrix_el                                                           */
/* ------------------------------------------------------------------ */

/* Step class reduction: DOWN folds into UP; EMPTY/FULL keep their identity.
 * Index by the raw step 1..4 (mirrors emol_expression.rb:173 STATE_CLASS). */
static const int STATE_CLASS[5] = { 0, EMPTY, UP, UP, FULL };

/* One row of the one-orbital matrix-element table (Table 4). The value is
 * mult * hexpr_root_val[root]; storing (mult,root) instead of a bare double keeps
 * this a compile-time `static const` (sqrt is not a constant expression) while
 * reproducing the previous branch-tree values bit-for-bit. */
typedef struct { int gc, hf, fc; double mult; unsigned char root; } matrix_el_row_t;

/* The 18 non-zero one-orbital matrix elements, class-reduced (gc,hf,fc).
 * Transcribed 1:1 from RubyRef/emol_expression.rb MATRIX_EL_TABLE (:175-195). */
static const matrix_el_row_t MATRIX_EL_TABLE[18] = {
    { EMPTY, OP_UNITY, EMPTY,  1, R1 },
    { EMPTY, OP_A,     UP,     1, R2 },
    { EMPTY, OP_AAS,   FULL,   1, R2 },
    { EMPTY, OP_ACS,   EMPTY, -1, R2 },
    { FULL,  OP_UNITY, FULL,   1, R1 },
    { FULL,  OP_C,     UP,    -1, R2 },
    { FULL,  OP_CAS,   FULL,  -1, R2 },
    { FULL,  OP_CCS,   EMPTY,  1, R2 },
    { FULL,  OP_CCA,   UP,     1, R1 },
    { FULL,  OP_CCAA,  FULL,   1, R1 },
    { UP,    OP_UNITY, UP,     1, R2 },
    { UP,    OP_A,     FULL,   1, R2 },
    { UP,    OP_C,     EMPTY,  1, R2 },
    { UP,    OP_CAS,   UP,    -1, R1 },
    { UP,    OP_CAT,   UP,     1, R3 },
    { UP,    OP_CAA,   FULL,  -1, R1 },
    { UP,    OP_ACS,   UP,    -1, R1 },
    { UP,    OP_ACT,   UP,     1, R3 },
};

/**
 * @brief Return the segment matrix-element factor for one walk-step pair under
 *        one operator: an 18-row table lookup (was a branch tree).
 * @param gf the "ket"/g-walk step code (@c EMPTY/UP/DOWN/FULL).
 * @param hf the operator code (@c OP_UNITY .. @c OP_ACT).
 * @param ff the "bra"/f-walk step code (@c EMPTY/UP/DOWN/FULL).
 * @return the one-orbital factor, or 0.0 when the combination is not in the table.
 *
 * Class-reduce gf/ff (DOWN->UP), scan @ref MATRIX_EL_TABLE for (gc,hf,fc), and
 * return @c mult*hexpr_root_val[root]. A miss returns 0.0 -- this reproduces the
 * old branch tree's (FULL,CAT,FULL)->0 silent zero AND its default-0 for any
 * invalid combination, byte-for-byte. Values are Table 4 (one-orbital matrix
 * elements <g~ h f>); see emol_expression.rb:197-212 for the same reduce+lookup.
 */
static double matrix_el(int gf, int hf, int ff)
{
    if (LK.phase == LK_SEARCH) LK.search_mel_calls++;  /* SEARCH must not reach here */
    int gc = STATE_CLASS[gf];
    int fc = STATE_CLASS[ff];
    for (int i = 0; i < 18; i++) {
        const matrix_el_row_t *r = &MATRIX_EL_TABLE[i];
        if (r->gc == gc && r->hf == hf && r->fc == fc)
            return r->mult * hexpr_root_val[r->root];
    }
    return 0.0;   /* (FULL,CAT,FULL) and all invalid combinations -> 0 */
}

/* ------------------------------------------------------------------ */
/* hmlkva                                                              */
/* ------------------------------------------------------------------ */

/**
 * @brief Evaluate one segment's full coupling factor: a phase times the
 *        @ref matrix_el value times a multiplicity root times a Wigner 9j.
 * @param gf,ff      the g-/f-walk step codes for this segment.
 * @param hf         the operator code.
 * @param psi,phi    the g-/f spin labels (2S+1 reps) at this node.
 * @param theta,theta_prev the operator's coupling label at this and the previous level.
 * @return the segment factor, or 0.0 if any reconstructed rep is non-positive.
 *
 * Role: called once per visited segment by @ref tree_search to extend the
 * running matrix-element product. It first reconstructs the previous-level spin
 * labels @c psi_prev / @c phi_prev by undoing the step (UP lowers, DOWN raises),
 * guards against non-positive reps, then assembles
 * @f$ \mathrm{phase}\cdot \mathrm{mel}\cdot\sqrt{\psi\,\phi_{prev}\,\theta}\cdot
 * \{9j\} @f$, where @c mel is @ref matrix_el and the 9j comes from @c wigner.h.
 * The phase is @f$ (-1)^{(f-1)(\theta_{prev}-1)} @f$ via the @c REPR_STATE /
 * @c REPR_OP rank tables.
 *
 * DOMAIN (per paper, Honda & Noro): This realizes the normalization-stripped
 * per-orbital factor L_j of Eq.18:
 *   L_j = (-1)^(n(Theta_{j-1}) n(f_j)) * sqrt([theta_j][Phi_{j-1}][Psi_j])
 *         * {9j-symbol} * <g~_j h_j f_j>.
 * The sqrt(...) here is the dimension factor with the PRE-step bra dimension
 * [Phi_{j-1}] (this is the L_j form actually implemented, not the symmetric
 * L'_j of Eq.13; see Sec. 2.5). The phase is the fermion sign
 * (-1)^(n(Theta)n(f)); the 9j-symbol couples the spin ranks of the bra, ket,
 * and operator chains at orbitals j-1 and j (Eq.13); and matrix_el() supplies
 * the one-orbital factor <g~ h f> (Table 4). The 9j value itself is evaluated in
 * wig/ (out of scope here); hmlkva assembles the Eq.18 product. See Eq.18,
 * Eq.13, Table 4, Sec. 2.5.
 */
static double hmlkva(int gf, int ff, int hf,
                     int psi, int phi, int theta, int theta_prev)
{
    int psi_prev, phi_prev;

    /* Reconstruct the previous-level spin labels by undoing this segment's step
     * (UP lowers by one, DOWN raises by one; EMPTY/FULL leave it unchanged). */
    if (gf == EMPTY || gf == FULL) psi_prev = psi;
    else if (gf == UP)             psi_prev = psi - 1;
    else                           psi_prev = psi + 1;   /* DOWN */

    if (ff == EMPTY || ff == FULL) phi_prev = phi;
    else if (ff == UP)             phi_prev = phi - 1;
    else                           phi_prev = phi + 1;   /* DOWN */

    /* Any non-positive rep means an impossible coupling -- contributes nothing. */
    if (psi <= 0 || psi_prev <= 0 || phi <= 0 || phi_prev <= 0 ||
        theta <= 0 || theta_prev <= 0)
        return 0.0;

    int h = REPR_OP[hf];
    int f = REPR_STATE[ff];
    int g = REPR_STATE[gf];
    int exp_val = (f - 1) * (theta_prev - 1);
    double phase = (exp_val & 1) ? -1.0 : 1.0;
    double mel   = matrix_el(gf, hf, ff);
    double mult  = sqrt((double)(psi * phi_prev * theta));
    if (LK.phase == LK_SEARCH) LK.search_wig_calls++;  /* SEARCH must not reach here */
    double w9j   = wigner_9j(h, theta_prev, theta,
                             f, phi_prev,   phi,
                             g, psi_prev,   psi);
    return phase * mel * mult * w9j;
}

/* ------------------------------------------------------------------ */
/* lk_value -- THE hot-line L(k) accessor (two explicit passes)        */
/* ------------------------------------------------------------------ */
/**
 * @brief Return L(k) for one segment: a table lookup, with an explicit
 *        compute-on-miss that is the fill path (GENERATE / OFF) or a counted
 *        safety net (SEARCH), never the intended search path.
 *
 *   HIT  -> pure lookup. In SEARCH this is THE path (no hmlkva/wigner/matrix_el).
 *   MISS -> GENERATE: expected -- warm is filling the table.
 *           OFF     : pair path -- fill-on-miss + cache, counted informationally.
 *           SEARCH  : a warm gap -- counted (miss + the wig/mel bumps inside
 *                     hmlkva); we assert this stays 0 on the whole oracle set.
 * hmlkva (=> matrix_el + wigner_9j) is thus reached ONLY on a miss.
 */
static double lk_value(int gf, int ff, int hf,
                       int psi, int phi, int theta, int theta_prev)
{
    uint32_t key = lk_pack(gf, ff, hf, psi, phi, theta, theta_prev);
    double v;
    if (lk_get(&LK, key, &v)) return v;                 /* HIT: pure lookup */
    if      (LK.phase == LK_SEARCH) LK.miss++;          /* counted safety net */
    else if (LK.phase == LK_OFF)    LK.off_fills++;     /* pair fill-on-miss */
    v = hmlkva(gf, ff, hf, psi, phi, theta, theta_prev);/* compute (fill / net) */
    lk_put(&LK, key, v);
    return v;
}

/* ------------------------------------------------------------------ */
/* Record hash for loop deduplication in tree_search                   */
/* Uses FNV-1a 64-bit on (bottom, top, max_w, min_w, ij0..3)          */
/* ------------------------------------------------------------------ */

#define RHASH_CAP 1024

typedef struct rnode { uint64_t h; struct rnode *next; } rnode_t;

typedef struct {
    rnode_t *buckets[RHASH_CAP];
    rnode_t  pool[256];   /* typical loops per call: small */
    int      pool_used;
} rhash_t;

/**
 * @brief Reset a record hash to empty (clear buckets, reset the node pool).
 * @param r record hash to initialize (caller-owned, e.g. a stack variable).
 * @note Does not free anything; pair a populated table with @ref rhash_reset to
 *       release any overflow heap nodes. Called once per start node in
 *       @ref calc_TensorOp_average.
 */
static void rhash_init(rhash_t *r)
{
    memset(r->buckets, 0, sizeof(r->buckets));
    r->pool_used = 0;
}

/**
 * @brief FNV-1a 64-bit hash over an int array.
 * @param arr values to hash.
 * @param n   number of values.
 * @return the 64-bit digest.
 * @note Used to key loops in @ref rhash_has_or_add for dedup.
 */
static uint64_t fnv64(const int *arr, int n)
{
    uint64_t h = 14695981039346656037ULL;
    for (int i = 0; i < n; i++) {
        h ^= (uint64_t)(unsigned int)arr[i];
        h *= 1099511628211ULL;
    }
    return h;
}

/**
 * @brief Test-and-insert a loop signature: report if this loop was already seen,
 *        otherwise record it.
 * @param r      record hash.
 * @param bottom,top loop endpoint node indices.
 * @param wg,wf  the g-/f path weights (folded order-independently into the key).
 * @param ijkl   orbital-index array of the loop.
 * @param nijkl  number of valid entries in @p ijkl.
 * @return 1 if the loop was already present (caller should skip it), 0 if it was
 *         new (now recorded).
 * @note The key is @c (bottom, top, max(wg,wf), min(wg,wf), ijkl...) so the two
 *       path weights are compared symmetrically. New nodes come from the inline
 *       pool until it is exhausted, then from @c malloc (freed by @ref rhash_reset).
 */
static int rhash_has_or_add(rhash_t *r, int bottom, int top,
                             int wg, int wf, const int *ijkl, int nijkl)
{
    int w_big = wg > wf ? wg : wf;
    int w_sml = wg > wf ? wf : wg;
    int key_arr[8];
    key_arr[0] = bottom; key_arr[1] = top;
    key_arr[2] = w_big;  key_arr[3] = w_sml;
    for (int i = 0; i < nijkl; i++) key_arr[4+i] = ijkl[i];
    uint64_t h = fnv64(key_arr, 4 + nijkl);
    uint32_t slot = (uint32_t)(h % RHASH_CAP);
    for (rnode_t *n = r->buckets[slot]; n; n = n->next)
        if (n->h == h) return 1;  /* found */
    /* insert */
    rnode_t *n;
    if (r->pool_used < 256)
        n = &r->pool[r->pool_used++];
    else
        n = malloc(sizeof(rnode_t));
    n->h = h;
    n->next = r->buckets[slot];
    r->buckets[slot] = n;
    return 0;  /* was new */
}

/* Cleanup rhash: free only heap-allocated nodes */
/**
 * @brief Clear a record hash, freeing only its heap-overflow nodes.
 * @param r          record hash to clear.
 * @param pool_limit number of inline-pool nodes in use (typically @c r->pool_used
 *                   captured before reset); nodes within @c [pool, pool+pool_limit)
 *                   are inline and must NOT be freed.
 * @note Frees only nodes outside the inline pool (those that came from @c malloc),
 *       then nulls all buckets. The pointer-range test is what distinguishes
 *       pooled from heap nodes. Counterpart to @ref rhash_init.
 */
static void rhash_reset(rhash_t *r, int pool_limit)
{
    for (int i = 0; i < RHASH_CAP; i++) {
        rnode_t *n = r->buckets[i];
        while (n) {
            rnode_t *next = n->next;
            if (n < r->pool || n >= r->pool + pool_limit)
                free(n);
            n = next;
        }
        r->buckets[i] = NULL;
    }
    r->pool_used = 0;
}

/* ------------------------------------------------------------------ */
/* Loop element vector                                                 */
/* ------------------------------------------------------------------ */

/**
 * @brief Append a loop element (by value copy) to a growable vector.
 * @param v target vector.
 * @param e element to copy in.
 * @note Stores a copy of @c *e (the vector does not alias caller storage). May
 *       @c realloc and move @c v->data, doubling capacity. Backing append used
 *       both per-node and for the merged result in @ref calc_TensorOp_average.
 */
static void le_push(le_vec_t *v, const loop_element_t *e)
{
    if (v->size == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 16;
        v->data = realloc(v->data, v->cap * sizeof(loop_element_t));
    }
    v->data[v->size++] = *e;
}

/**
 * @brief qsort/stable-sort comparator ordering ::loop_element_t by @c top.
 * @param a,b pointers to ::loop_element_t.
 * @return a->top - b->top.
 * @note Used with @c hexpr_stable_sort so equal-top elements keep walk order.
 */
static int le_cmp_top(const void *a, const void *b)
{
    return ((const loop_element_t *)a)->top - ((const loop_element_t *)b)->top;
}

/* ------------------------------------------------------------------ */
/* tree_search (recursive)                                             */
/* ------------------------------------------------------------------ */

/* State passed into (and mutated by) the recursion; cleaned up after
 * each segment via push/pop pattern.  Stack lives on caller's frame. */

/**
 * @brief Recursively extend a double walk (g and f) down the DRT by one segment,
 *        emitting a ::loop_element_t whenever a loop closes.
 * @param drt     the DRT being walked (read via the @c hexpr_drt_* accessors).
 * @param bottom  the start node index shared by every loop on this descent.
 * @param gf,ff   the g-/f step codes (1..4) taken into this level.
 * @param k       the current orbital level (the recursion descends to @c k+1).
 * @param ng,nf   per-level g-/f node-index stacks, valid @c [0..norb].
 * @param wg,wf   per-level g-/f cumulative path weights.
 * @param iop     index of the current operator within @p tp.
 * @param tp      the tensor product driving the walk.
 * @param hf      the current operator code (selects @c NUM_CA and the segment factor).
 * @param theta,theta_prev the operator coupling labels at this/previous level.
 * @param mk      per-level running product of segment factors.
 * @param ijkl    orbital-index stack; @param nijkl its current depth (both mutated).
 * @param origin  tag copied into every emitted loop element (set by the caller).
 * @param expr    output loop-element vector for this start node.
 * @param record  dedup hash so identical loops are emitted once.
 * @param norb    orbital count (bounds the descent).
 *
 * Role: the core DRT walk. Each call consumes one segment: it follows the g and
 * f upper arcs for steps @c gf/ff, accumulates path weights and the
 * @ref hmlkva segment factor into the running product @c mk, then either closes
 * the loop or recurses deeper. The per-level state arrays act as an explicit
 * stack indexed by orbital level, and the @c ijkl push at entry / pop at exit
 * keeps the orbital-index stack consistent across branches.
 *
 * Control flow: dead branches (a missing arc or weight) return immediately;
 * below the product threshold @c THRESH_ZERO the subtree is abandoned. A loop
 * closes only when the LAST operator has been placed (@c iop==lenop-1) and the
 * two walks rejoin (@c ng==nf) -- then a deduped loop element is pushed.
 * Otherwise the walk continues two ways: inserting @c OP_UNITY filler segments
 * (guarded so enough orbitals remain for the operators still to place) and/or
 * advancing to the next operator's segment-walk set.
 *
 * DOMAIN (per paper, Honda & Noro): The recursion walks the DRT one orbital at a
 * time, calling hmlkva() at each level to get the per-orbital factor L_j and
 * accumulating their product -- this is the product form prod_j L_j of Eq.16
 * (equivalently Eq.12, the reduction of the N-electron matrix element to
 * per-orbital factors). A "loop" is the segment of the walk where the bra (Psi)
 * and ket (Phi) paths diverge from a common DRT node and rejoin at another; it
 * is "closed" when the two walks coincide at both ends, so that
 * [Phi_{a-1}] = [Phi_b] and the per-orbital normalization telescopes to unity
 * over the loop (Eq.19, Eq.20). Outside a loop the bra and ket walks are
 * identical and contribute the trivial factor. See Sec. 2.4-2.5, Eq.16.
 */
static void tree_search(const hexpr_drt_t *drt, int bottom,
                        int gf, int ff,       /* current walk pair */
                        int k,                /* current orb level */
                        int *ng, int *nf,     /* node arrays [0..norb] */
                        int *wg, int *wf,     /* path weight arrays */
                        int iop,              /* current operator index */
                        const tensor_product_t *tp,
                        int hf,               /* current operator code */
                        int theta, int theta_prev,
                        double *mk,           /* matrix element products */
                        int *ijkl, int *nijkl, /* orbital index stack */
                        int origin,
                        le_vec_t *expr, rhash_t *record,
                        int norb)
{
    /* Follow the g and f upper arcs for this step pair; a missing arc or weight
     * (negative) means this branch is dead -- stop before mutating any state. */
    int next_ng = hexpr_drt_upper_arc(drt, ng[k], gf - 1);
    int next_nf = hexpr_drt_upper_arc(drt, nf[k], ff - 1);
    int aw_g    = hexpr_drt_arc_weight(drt, ng[k], gf - 1);
    int aw_f    = hexpr_drt_arc_weight(drt, nf[k], ff - 1);

    if (next_ng < 0 || next_nf < 0 || aw_g < 0 || aw_f < 0)
        return;

    /* Push level k+1: record both walk nodes and the cumulative path weights. */
    int k1 = k + 1;
    ng[k1] = next_ng;
    nf[k1] = next_nf;
    wg[k1] = wg[k] + aw_g;
    wf[k1] = wf[k] + aw_f;

    /* This operator contributes num_ca orbital indices (its creation+annihilation
     * count); push them onto the shared ijkl stack at this level. */
    int num_ca = NUM_CA[hf];
    for (int i = 0; i < num_ca; i++)
        ijkl[(*nijkl)++] = k1;

    /* Extend the running matrix-element product by this segment's factor.
     * Two-phase L(k): pure table lookup on the search pass; the table was filled
     * by warm_lk_table in the generate pass (see expr.c). hmlkva runs only on a
     * table miss (fill in generate/off; counted safety net in search). */
    int psi = hexpr_drt_spin(drt, ng[k1]);
    int phi = hexpr_drt_spin(drt, nf[k1]);
    mk[k1] = mk[k] * lk_value(gf, ff, hf, psi, phi, theta, theta_prev);

    /* Below threshold the subtree cannot contribute -- prune it (still pops). */
    if (fabs(mk[k1]) >= THRESH_ZERO) {
        if (iop == tp->lenop - 1 && ng[k1] == nf[k1]) {
            /* loop closed */
            /* All operators placed and the walks rejoined: emit one loop element,
             * unless this exact loop signature was already recorded. */
            if (!rhash_has_or_add(record, bottom, ng[k1],
                                  wg[k1], wf[k1], ijkl, *nijkl)) {
                loop_element_t le;
                le.bottom   = bottom;
                le.top      = ng[k1];
                le.weight_g = wg[k1];
                le.weight_f = wf[k1];
                le.num_ijkl = *nijkl;
                for (int i = 0; i < *nijkl; i++) le.ijkl[i] = ijkl[i];
                le.value    = mk[k1];
                le.seq_num  = 0;
                le.origin   = origin;
                le_push(expr, &le);
            }
        } else {
            /* UNITY segments between operators */
            /* Continue with filler identity segments, but only while enough
             * orbitals remain below for the operators still to be placed. */
            if (norb - k1 > tp->lenop - iop - 1) {
                const wset_t *ws = &SEG_WALK[OP_UNITY];
                for (int w = 0; w < ws->n; w++)
                    tree_search(drt, bottom,
                                ws->w[w][0], ws->w[w][1],
                                k1, ng, nf, wg, wf,
                                iop, tp, OP_UNITY, theta, theta, mk,
                                ijkl, nijkl, origin, expr, record, norb);
            }
            /* advance to next operator */
            /* Place the next operator: its code/theta drive a new segment-walk
             * set, and the current theta becomes theta_prev for that level. */
            if (tp->lenop - 1 > iop) {
                int next_iop   = iop + 1;
                int next_hf    = tp->op[next_iop];
                int next_theta = tp->theta[next_iop];
                const wset_t *ws = &SEG_WALK[next_hf];
                for (int w = 0; w < ws->n; w++)
                    tree_search(drt, bottom,
                                ws->w[w][0], ws->w[w][1],
                                k1, ng, nf, wg, wf,
                                next_iop, tp, next_hf, next_theta, theta,
                                mk, ijkl, nijkl, origin, expr, record, norb);
            }
        }
    }

    /* pop */
    /* Restore the ijkl stack depth so sibling branches start clean (the per-level
     * ng/nf/wg/wf/mk slots are simply overwritten by the next branch). */
    *nijkl -= num_ca;
}

/* ------------------------------------------------------------------ */
/* calc_TensorOp_average                                               */
/* ------------------------------------------------------------------ */

/* Returns flat array of loop elements (caller frees le_vec_t.data).
 * Result is sorted by .top, duplicates already removed by record hash. */
/**
 * @brief Generate all loop elements for a tensor product over a DRT.
 * @param drt    the DRT to walk (from @c drt.c).
 * @param tp     the tensor product (operators + theta labels) to evaluate.
 * @param origin a caller tag stamped onto every produced loop element.
 * @return a flat ::le_vec_t of loop elements; the CALLER owns and must free
 *         @c .data. Within each start node the elements are sorted by @c top and
 *         de-duplicated; @c .value carries the loop's matrix-element product.
 *
 * Role: the public entry point of this unit and the function @c expr.c's
 * @c caseNN routines invoke (via the @c CALC macro). It drives @ref tree_search
 * from every admissible bottom node and concatenates the per-node results.
 *
 * Structure: allocate the per-level walk-state arrays once (sized @c norb+1),
 * then for each start node that still leaves room for all @c lenop operators
 * (the @c norb-orb_i cutoff), seed level-@c orb_i state (both walks at that node,
 * zero weights, unit product), launch @ref tree_search over the first operator's
 * segment-walk set, stable-sort the node's loops by @c top, and append them to
 * the running result. The record hash is re-initialized per node and its heap
 * overflow reclaimed by @ref rhash_reset before moving on.
 *
 * DOMAIN (per paper, Honda & Noro): "average" here is a historical name (from
 * the Ruby prototype), not an arithmetic mean. The routine walks the DRT from
 * each admissible start node (via tree_search), collecting the per-orbital
 * product prod_j L_j for every loop that completes, and concatenates the
 * results (sorted by .top). The normalization of the generator-form CSFs is not
 * divided out here: the genealogical CSFs are un-normalized (norm sqrt([Phi]),
 * Eq.15), but the L_j accumulated by hmlkva already carry the normalization
 * split, using the pre-step bra dimension [Phi_{j-1}] (the implemented form of
 * Eq.16/Eq.18, not the symmetric L'_j of Eq.13). The concatenated set of
 * per-(bra,ket) coupling contributions is what the CALC callers in expr.c /
 * csf_pair_core.c consume. See Sec. 2.5, Eq.15, Eq.16, Eq.18.
 */
le_vec_t calc_TensorOp_average(const hexpr_drt_t *drt,
                                const tensor_product_t *tp, int origin)
{
    int nnodes = hexpr_drt_nnodes(drt);
    int norb   = hexpr_drt_norb(drt);
    int lenop  = tp->lenop;

    /* per-level state arrays, indexed 0..norb */
    int *ng = malloc((norb+1) * sizeof(int));
    int *nf = malloc((norb+1) * sizeof(int));
    int *wg = malloc((norb+1) * sizeof(int));
    int *wf = malloc((norb+1) * sizeof(int));
    double *mk = malloc((norb+1) * sizeof(double));

    int  ijkl[4];
    int  nijkl;

    le_vec_t all_expr = {NULL, 0, 0};
    rhash_t  record;

    int iop = 0;
    int hf0 = tp->op[iop];

    for (int inode = 0; inode < nnodes; inode++) {
        /* Stop once too few orbitals remain above this node to place all
         * lenop operators -- no loop could complete from here. */
        int orb_i = hexpr_drt_orb(drt, inode);
        if (norb - orb_i <= lenop - 1) break; /* no more room */

        le_vec_t expr = {NULL, 0, 0};
        rhash_init(&record);

        /* Seed both walks at this start node: same node, zero path weights,
         * unit product, empty orbital-index stack. */
        int k = orb_i;
        ng[k] = inode; nf[k] = inode;
        wg[k] = 0;     wf[k] = 0;
        mk[k] = 1.0;
        nijkl  = 0;

        int theta      = tp->theta[iop];
        int theta_prev = tp->theta_prev[iop]; /* = 1 */

        const wset_t *ws = &SEG_WALK[hf0];
        for (int w = 0; w < ws->n; w++)
            tree_search(drt, inode,
                        ws->w[w][0], ws->w[w][1],
                        k, ng, nf, wg, wf,
                        iop, tp, hf0, theta, theta_prev,
                        mk, ijkl, &nijkl, origin, &expr, &record, norb);

        /* sort by .top (ascending) -- stable sort, matches Ruby x.top-y.top */
        hexpr_stable_sort(expr.data, expr.size, sizeof(loop_element_t), le_cmp_top);

        /* append to all_expr */
        for (int i = 0; i < expr.size; i++)
            le_push(&all_expr, &expr.data[i]);

        /* Release this node's scratch; rhash_reset frees only heap overflow. */
        free(expr.data);
        rhash_reset(&record, record.pool_used);
    }

    free(ng); free(nf); free(wg); free(wf); free(mk);
    return all_expr;
}
