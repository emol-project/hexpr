/* loopgen/brooks_table.h -- declarative Brooks case catalog (data only).
 *
 * THE single source of truth for Brooks' case catalog, driving BOTH the
 * whole-matrix path (expr.c brooks_all) and the per-pair path
 * (csf_pair_core.c pair_eval_internal). Reading a row tells you the whole
 * case; interpreting it needs NO per-case branching and NO function pointers.
 *
 * CROSS-CHECK (field-for-field, verified in Phase A-1): every row agrees with
 *   RubyRef/emol_expression.rb BROOKS_CASES (:664-819) [source of truth],
 *   expr.c case01..case14, and csf_pair_core.c case01_pair..case14_pair --
 *   id, origin, enable_whole, enable_pair, each op's chain+theta, each term's
 *   op_indices+coef+perm.
 *
 * NOTE (case9): row present with enable_whole=false AND enable_pair=false, so
 *   it NEVER emits -- both drivers skip it via enable_*. Kept for a 1:1 mapping
 *   with the Ruby table; the enabled-row sequence still equals brooks_all /
 *   pair_eval_internal exactly.
 *
 * NOTE (coef encoding): a coefficient is stored as {mult, root} and resolved to
 *   mult * hexpr_root_val[root] (root: R1=>*1, R2=>*sqrt2, R3=>*sqrt3). This keeps
 *   the table a compile-time `static const` (sqrt is not a constant expression)
 *   while reproducing the Ruby doubles bit-for-bit (mult is an exact rational:
 *   +-1, +2, +0.5). hexpr_root_val is the single shared roots array (expr_internal.h).
 */
#ifndef HEXPR_BROOKS_TABLE_H
#define HEXPR_BROOKS_TABLE_H

#include "expr_internal.h"   /* EMPTY..FULL, OP_*, R1/R2/R3, hexpr_root_val */

/* Actual maxima across the 17 cases: nops<=2, lenop<=4, nterms<=2,
 * n_op_indices<=2, perm_len in {2,4}. op_indices/coef sized 4 per field spec. */
#define BC_MAX_OPS    2
#define BC_MAX_TERMS  2
#define BC_MAX_CHAIN  4
#define BC_MAX_OPIDX  4
#define BC_MAX_PERM   4

typedef struct {              /* one coefficient = mult * sqrt(root) */
    double        mult;       /* exact rational: +-1, +2, +0.5 */
    unsigned char root;       /* R1 / R2 / R3 */
} brooks_coef_t;

typedef struct {              /* one tensor product (chain of operators) */
    int chain[BC_MAX_CHAIN];  /* OP_* codes; unused slots = 0 */
    int theta[BC_MAX_CHAIN];  /* coupling labels; unused slots = 0 */
    int lenop;                /* active length of chain/theta (1..4) */
} brooks_op_t;

typedef struct {              /* one add_expr / add_expr_pair call */
    int           op_indices[BC_MAX_OPIDX]; /* which ops feed this term */
    brooks_coef_t coef[BC_MAX_OPIDX];       /* coef[i] scales op_indices[i] */
    int           n_op_indices;             /* active length (1 or 2) */
    int           perm[BC_MAX_PERM];        /* index map handed to add_expr */
    int           perm_len;                 /* 2 (1e) or 4 (2e) */
} brooks_term_t;

typedef struct {
    const char   *id;            /* label only; NEVER dispatched on */
    int           origin;        /* 3rd arg to calc_TensorOp_average (CALC) */
    int           enable_whole;  /* included by the whole-matrix driver */
    int           enable_pair;   /* included by the per-pair driver */
    brooks_op_t   ops[BC_MAX_OPS];
    int           nops;
    brooks_term_t terms[BC_MAX_TERMS];
    int           nterms;
} brooks_case_t;

/* Row order == Ruby BROOKS_CASES order: case01,02,03, case1..14 (case9 disabled).
 * Enabled-row sequence reproduces the current emission order exactly:
 *   whole (skip !enable_whole): case01,02,1..8,10..14   == brooks_all
 *   pair  (skip !enable_pair):  case01,02,03,1..8,10..14 == pair_eval_internal */
static const brooks_case_t BROOKS_CASES[] = {
  /* --- one-electron cases (Table 3) --- */
  { "case01", 0, 1, 1,
    { {{OP_CAS,0,0,0},{1,0,0,0},1} }, 1,
    { { {0},{{-1,R2}},              1, {0,1},   2 } }, 1 },

  { "case02", 0, 1, 1,
    { {{OP_C,OP_A,0,0},{2,1,0,0},2} }, 1,
    { { {0},{{-1,R2}},              1, {0,1},   2 } }, 1 },

  { "case03", 0, 0, 1,                            /* enable_whole=false */
    { {{OP_A,OP_C,0,0},{2,1,0,0},2} }, 1,
    { { {0},{{-1,R2}},              1, {0,1},   2 } }, 1 },

  /* --- two-electron cases (Table 1/2) --- */
  { "case1", 1, 1, 1,
    { {{OP_A,OP_C,OP_C,OP_A},{2,1,2,1},4},
      {{OP_A,OP_C,OP_C,OP_A},{2,3,2,1},4} }, 2,
    { { {0,1},{{-1,R1},{-1,R3}},   2, {0,2,1,3},4 },   /* ikjl */
      { {0},  {{2,R1}},            1, {0,1,2,3},4 } }, 2 }, /* ijkl */

  { "case2", 2, 1, 1,
    { {{OP_C,OP_A,OP_C,OP_A},{2,1,2,1},4},
      {{OP_C,OP_A,OP_C,OP_A},{2,3,2,1},4} }, 2,
    { { {0,1},{{-1,R1},{1,R3}},    2, {0,3,1,2},4 },   /* iljk */
      { {0},  {{2,R1}},            1, {0,1,2,3},4 } }, 2 }, /* ijkl */

  { "case3", 3, 1, 1,
    { {{OP_C,OP_C,OP_A,OP_A},{2,1,2,1},4},
      {{OP_C,OP_C,OP_A,OP_A},{2,3,2,1},4} }, 2,
    { { {0,1},{{-1,R1},{-1,R3}},   2, {0,2,1,3},4 },   /* ikjl */
      { {0,1},{{-1,R1},{1,R3}},    2, {0,3,1,2},4 } }, 2 }, /* iljk */

  { "case4", 4, 1, 1,
    { {{OP_A,OP_CCS,OP_A,0},{2,2,1,0},3} }, 1,
    { { {0},{{1,R1}},              1, {0,1,2,3},4 } }, 1 },

  { "case5", 5, 1, 1,
    { {{OP_C,OP_CAS,OP_A,0},{2,2,1,0},3},
      {{OP_C,OP_CAT,OP_A,0},{2,2,1,0},3} }, 2,
    { { {0,1},{{-1,R1},{1,R3}},    2, {0,1,2,3},4 },   /* ijjl */
      { {0},  {{2,R1}},            1, {0,3,1,2},4 } }, 2 }, /* jjil */

  { "case6", 6, 1, 1,
    { {{OP_CCS,OP_A,OP_A,0},{1,2,1,0},3} }, 1,
    { { {0},{{1,R1}},              1, {0,3,1,2},4 } }, 1 },

  { "case7", 7, 1, 1,
    { {{OP_CAS,OP_C,OP_A,0},{1,2,1,0},3},
      {{OP_CAT,OP_C,OP_A,0},{3,2,1,0},3} }, 2,
    { { {0,1},{{-1,R1},{-1,R3}},   2, {0,2,1,3},4 },   /* ikil */
      { {0},  {{2,R1}},            1, {0,1,2,3},4 } }, 2 }, /* iikl */

  { "case8", 8, 1, 1,
    { {{OP_C,OP_C,OP_AAS,0},{2,1,1,0},3} }, 1,
    { { {0},{{-1,R1}},             1, {0,3,1,2},4 } }, 1 },

  { "case9", 9, 0, 0,               /* DISABLED in both paths; never emits.  */
    { {{OP_A,OP_C,OP_CAS,0},{2,1,1,0},3},  /* present only to keep the table   */
      {{OP_A,OP_C,OP_CAT,0},{2,3,1,0},3} }, 2, /* 1:1 with Ruby BROOKS_CASES.   */
    { { {0,1},{{-1,R1},{1,R3}},    2, {0,2,1,3},4 } }, 1 }, /* iljl */

  { "case10", 10, 1, 1,
    { {{OP_C,OP_A,OP_CAS,0},{2,1,1,0},3},
      {{OP_C,OP_A,OP_CAT,0},{2,3,1,0},3} }, 2,
    { { {0,1},{{-1,R1},{-1,R3}},   2, {0,2,1,3},4 },   /* iljl */
      { {0},  {{2,R1}},            1, {0,1,2,3},4 } }, 2 }, /* ijll */

  { "case11", 11, 1, 1,
    { {{OP_C,OP_CAA,0,0},{2,1,0,0},2},
      {{OP_CCA,OP_A,0,0},{2,1,0,0},2} }, 2,
    { { {0},{{2,R1}},              1, {0,1,2,3},4 },   /* illl (op0 only) */
      { {1},{{2,R1}},              1, {0,1,2,3},4 } }, 2 }, /* iiil (op1 only) */

  { "case12", 12, 1, 1,
    { {{OP_CCS,OP_AAS,0,0},{1,1,0,0},2} }, 1,
    { { {0},{{0.5,R1}},            1, {0,2,1,3},4 } }, 1 },

  { "case13", 13, 1, 1,
    { {{OP_CAS,OP_CAS,0,0},{1,1,0,0},2},
      {{OP_CAT,OP_CAT,0,0},{3,1,0,0},2} }, 2,
    { { {0,1},{{-1,R1},{1,R3}},    2, {0,2,1,3},4 },   /* ilil */
      { {0},  {{2,R1}},            1, {0,1,2,3},4 } }, 2 }, /* iill */

  { "case14", 14, 1, 1,
    { {{OP_CCAA,0,0,0},{1,0,0,0},1} }, 1,
    { { {0},{{1,R1}},              1, {0,1,2,3},4 } }, 1 },
};
#define N_BROOKS_CASES ((int)(sizeof(BROOKS_CASES)/sizeof(BROOKS_CASES[0])))
/* = 17 rows: case01,02,03, case1..14 (case9 present but disabled) -- 1:1 Ruby. */

/* Materialize a term's coefficients as mult*hexpr_root_val[root] into a caller
 * buffer (length >= tm->n_op_indices). Shared by the whole and pair drivers so
 * both resolve coefficients through the one roots array. */
static inline void bc_fill_coef(const brooks_term_t *tm, double *out)
{
    for (int i = 0; i < tm->n_op_indices; i++)
        out[i] = tm->coef[i].mult * hexpr_root_val[tm->coef[i].root];
}

#endif  /* HEXPR_BROOKS_TABLE_H */
