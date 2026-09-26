/* loopgen/expr_internal.h -- shared internal types for the expression generator */
#ifndef HEXPR_EXPR_INTERNAL_H_
#define HEXPR_EXPR_INTERNAL_H_

#include <stdint.h>
#include <stdlib.h>
#include "drt.h"
#include "sort.h"

/* ---- Walk step + operator codes (from emol_const_ops.rb) ----
 * Single shared definition; expr.c / expr_tree.c / csf_pair_core.c all use these. */
enum { EMPTY = 1, UP = 2, DOWN = 3, FULL = 4 };
enum {
    OP_UNITY = 1, OP_A = 2, OP_C = 3, OP_AAS = 4, OP_CAS = 5, OP_CAT = 6,
    OP_CCS = 7, OP_CAA = 8, OP_CCA = 9, OP_CCAA = 10, OP_ACS = 11, OP_ACT = 12
};

/* ---- sqrt-root multiplier codes for the {mult,root} coefficient encoding ----
 * A coefficient is stored as (mult, root) and resolved to mult * hexpr_root_val[root].
 * R1 => *1, R2 => *sqrt2, R3 => *sqrt3. */
enum { R1 = 1, R2 = 2, R3 = 3 };

/* Single shared roots array -- THE one sqrt(2)/sqrt(3) path in the tree.
 * hexpr_root_val = { unused, 1.0, sqrt(2.0), sqrt(3.0) }. hexpr_root_val[R1]=1.0
 * is exact from the static initializer; [R2]/[R3] are seeded by hexpr_roots_init().
 * Seeded at each driver's init (bc_init_constants / init_constants) before any walk. */
extern double hexpr_root_val[4];
void hexpr_roots_init(void);

/* ---- two-phase L(k) table: control + instrumentation (defined in expr_tree.c) ----
 * The whole path runs generate(warm) -> search; the pair path stays in OFF and
 * fills-on-miss. See expr_tree.c / expr.c for the two explicit passes. */
void hexpr_lk_ensure_init(void);      /* one-time clear of the shared table */
void hexpr_lk_begin_generation(void); /* ensure_init + phase = GENERATE (fill) */
void hexpr_lk_begin_search(void);     /* phase = SEARCH, reset search counters */
void hexpr_lk_finish(void);           /* phase = OFF */
int  hexpr_lk_phase(void);            /* current phase (LK_OFF/GENERATE/SEARCH) */
long hexpr_lk_search_misses(void);    /* SEARCH-phase table misses  (assert 0) */
long hexpr_lk_search_wig_calls(void); /* SEARCH-phase wigner_9j calls (assert 0) */
long hexpr_lk_search_mel_calls(void); /* SEARCH-phase matrix_el calls (assert 0) */
long hexpr_lk_off_fills(void);        /* OFF-phase (pair) fill-on-miss count */
long hexpr_lk_table_size(void);       /* live keys in the shared table */

/* ---- Tensor product (from Ruby TensorProduct) ---- */
typedef struct {
    int op[4];
    int theta[4];
    int theta_prev[4];
    int lenop;
} tensor_product_t;

static inline void tp_init(tensor_product_t *tp,
                            const int *op, const int *theta, int lenop)
{
    tp->lenop = lenop;
    for (int i = 0; i < lenop; i++) {
        tp->op[i]         = op[i];
        tp->theta[i]      = theta[i];
        tp->theta_prev[i] = (i == 0) ? 1 : theta[i-1];
    }
}

/* ---- Loop element (from Ruby LoopElement) ---- */
typedef struct {
    int    bottom, top;
    int    weight_g, weight_f;
    int    ijkl[4];
    int    num_ijkl;
    double value;
    int    seq_num;
    int    origin;
} loop_element_t;

/* ---- Dynamic array of loop_element_t ---- */
typedef struct {
    loop_element_t *data;
    int             size;
    int             cap;
} le_vec_t;

/* ---- Compressed expression term ---- */
typedef struct {
    int    ij_key;    /* triangle(lw[0]+uw[0]+wg, lw[0]+uw[0]+wf) */
    int    seq_num;   /* pqrs orbital index */
    int    bottom, top;
    int    weight_g, weight_f;
    double value;
} comp_term_t;

/* ---- Dynamic array of comp_term_t ---- */
typedef struct {
    comp_term_t *data;
    int          size;
    int          cap;
} ct_vec_t;

static inline void ct_push(ct_vec_t *v, const comp_term_t *t)
{
    if (v->size == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 64;
        v->data = realloc(v->data, v->cap * sizeof(comp_term_t));
    }
    v->data[v->size++] = *t;
}

/* ---- Functions implemented in expr_tree.c ---- */
le_vec_t calc_TensorOp_average(const hexpr_drt_t *drt,
                                const tensor_product_t *tp, int origin);

/* ---- Final expanded term (for output) ---- */
typedef struct {
    /* 0-based packed triangle index of the (bra,ket) pair. The binary APIs emit
     * this raw (0-based); the text dump alone adds +1 to match the Ruby corpus
     * reference/<case>.expr.txt (the intentional binary-0 / dump-1 duality). */
    int    ij;
    int    pqrs;  /* seq_num: 0-based packed integral index, emitted raw everywhere */
    double coef;
} out_term_t;

/* ---- Dynamic array of out_term_t ---- */
typedef struct {
    out_term_t *data;
    int         size;
    int         cap;
} ot_vec_t;

static inline void ot_push(ot_vec_t *v, const out_term_t *t)
{
    if (v->size == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 256;
        v->data = realloc(v->data, v->cap * sizeof(out_term_t));
    }
    v->data[v->size++] = *t;
}

#endif /* HEXPR_EXPR_INTERNAL_H_ */
