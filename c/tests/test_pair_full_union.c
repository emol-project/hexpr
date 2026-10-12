/* test_pair_full_union.c -- the v1.1 complete-decomposition calls.
 *
 * For every pair (I, J) of several DRTs (SOCI, FOCI, full CAS, singlet and
 * triplet, and CSFSET sub-DRTs) this checks, against the v1.0 calls:
 *
 *   1. hexpr_pair_eval_full(I,J) and (J,I) are the same term set with the
 *      same coefficients (order independence);
 *   2. that set equals the union of hexpr_pair_eval(I,J) and (J,I), with
 *      agreeing coefficients (docs/PAIR_API_SEMANTICS.md, G4/G5);
 *   3. for I == J it equals hexpr_pair_eval(I,I);
 *   4. no (ij, pqrs) key repeats within one output;
 *   5. hexpr_pair_count_full matches, and HEXPR_EXPR_ONE is the one-electron
 *      part of HEXPR_EXPR_FULL;
 *   6. hexpr_block_eval_full / _count_full over a bra x ket block emit, per
 *      ordered pair, exactly the hexpr_pair_eval_full terms.
 *
 * Large spaces are sampled (a fixed pseudo-random set of pairs). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "../include/hexpr.h"

typedef struct { int32_t ij, pqrs; double coef; } term_t;
typedef struct { term_t *t; int n; } tset_t;

static int failures = 0;
#define FAIL(...) do { if (failures < 20) { printf("FAIL: "); printf(__VA_ARGS__); printf("\n"); } failures++; } while (0)

static int term_cmp(const void *a, const void *b)
{
    const term_t *x = a, *y = b;
    if (x->ij != y->ij) return x->ij < y->ij ? -1 : 1;
    if (x->pqrs != y->pqrs) return x->pqrs < y->pqrs ? -1 : 1;
    return 0;
}

/* one call -> sorted term set */
static tset_t get(const hexpr_drt_t *drt, int I, int J, hexpr_expr_which_t w, int full)
{
    tset_t s = {NULL, 0};
    int mx = 0;
    hexpr_pair_max_terms(drt, w, &mx);
    int32_t *ij = malloc((size_t)(mx + 1) * sizeof *ij);
    int32_t *pq = malloc((size_t)(mx + 1) * sizeof *pq);
    double *cf = malloc((size_t)(mx + 1) * sizeof *cf);
    hexpr_status_t st = full ? hexpr_pair_eval_full(drt, I, J, w, &s.n, ij, pq, cf)
                             : hexpr_pair_eval(drt, I, J, w, &s.n, ij, pq, cf);
    if (st != HEXPR_OK) FAIL("pair_eval%s(%d,%d) status %d", full ? "_full" : "", I, J, st);
    s.t = malloc((size_t)(s.n + 1) * sizeof *s.t);
    for (int k = 0; k < s.n; k++) s.t[k] = (term_t){ij[k], pq[k], cf[k]};
    qsort(s.t, (size_t)s.n, sizeof *s.t, term_cmp);
    free(ij); free(pq); free(cf);
    return s;
}

static int has_dup(const tset_t *s)
{
    for (int k = 1; k < s->n; k++)
        if (term_cmp(&s->t[k - 1], &s->t[k]) == 0) return 1;
    return 0;
}

/* sorted merge of a and b; returns 0 if a shared key carries different coefs */
static tset_t merge(const tset_t *a, const tset_t *b, int *agree)
{
    tset_t m = {malloc((size_t)(a->n + b->n + 1) * sizeof(term_t)), 0};
    int i = 0, j = 0;
    *agree = 1;
    while (i < a->n || j < b->n) {
        int c = i >= a->n ? 1 : j >= b->n ? -1 : term_cmp(&a->t[i], &b->t[j]);
        if (c < 0) m.t[m.n++] = a->t[i++];
        else if (c > 0) m.t[m.n++] = b->t[j++];
        else {
            if (fabs(a->t[i].coef - b->t[j].coef) > 1e-12) *agree = 0;
            m.t[m.n++] = a->t[i++]; j++;
        }
    }
    return m;
}

static int same(const tset_t *a, const tset_t *b)
{
    if (a->n != b->n) return 0;
    for (int k = 0; k < a->n; k++)
        if (term_cmp(&a->t[k], &b->t[k]) != 0 || fabs(a->t[k].coef - b->t[k].coef) > 1e-12)
            return 0;
    return 1;
}

static long n_strict = 0;   /* pairs where one v1.0 direction was incomplete */

static void check_pair(const hexpr_drt_t *drt, const char *name, int I, int J)
{
    int n_oel = hexpr_drt_norb(drt);
    n_oel = n_oel * (n_oel + 1) / 2;
    tset_t A = get(drt, I, J, HEXPR_EXPR_FULL, 0), B = get(drt, J, I, HEXPR_EXPR_FULL, 0);
    tset_t F1 = get(drt, I, J, HEXPR_EXPR_FULL, 1), F2 = get(drt, J, I, HEXPR_EXPR_FULL, 1);
    if (has_dup(&F1) || has_dup(&F2)) FAIL("%s (%d,%d): repeated key in _full output", name, I, J);
    if (!same(&F1, &F2)) FAIL("%s (%d,%d): _full depends on argument order (%d vs %d terms)",
                              name, I, J, F1.n, F2.n);
    int agree = 1;
    tset_t U = merge(&A, &B, &agree);
    if (!agree) FAIL("%s (%d,%d): v1.0 directions disagree on a shared key", name, I, J);
    if (!same(&U, &F1)) FAIL("%s (%d,%d): _full (%d terms) != union of directions (%d terms)",
                             name, I, J, F1.n, U.n);
    if (I == J && !same(&A, &F1)) FAIL("%s (%d,%d): diagonal _full != pair_eval", name, I, J);
    if (I != J && A.n != B.n) n_strict++;

    int cnt = -1;
    hexpr_pair_count_full(drt, I, J, HEXPR_EXPR_FULL, &cnt);
    if (cnt != F1.n) FAIL("%s (%d,%d): pair_count_full %d != %d", name, I, J, cnt, F1.n);

    tset_t O = get(drt, I, J, HEXPR_EXPR_ONE, 1);
    int k1 = 0, ok = 1;
    for (int k = 0; k < F1.n; k++) {
        if (F1.t[k].pqrs >= n_oel) continue;
        if (k1 >= O.n || term_cmp(&F1.t[k], &O.t[k1]) != 0 || F1.t[k].coef != O.t[k1].coef) ok = 0;
        k1++;
    }
    if (!ok || k1 != O.n) FAIL("%s (%d,%d): ONE is not the one-electron part of FULL", name, I, J);
    hexpr_pair_count_full(drt, I, J, HEXPR_EXPR_ONE, &cnt);
    if (cnt != O.n) FAIL("%s (%d,%d): pair_count_full(ONE) %d != %d", name, I, J, cnt, O.n);

    free(A.t); free(B.t); free(F1.t); free(F2.t); free(U.t); free(O.t);
}

/* block_eval_full over bra x ket == per-pair _full, tagged with positions */
static void check_block(const hexpr_drt_t *drt, const char *name,
                        const int *bra, int nb, const int *ket, int nk)
{
    int mx = 0, cnt = -1, n = 0;
    hexpr_block_max_terms(drt, nb, nk, HEXPR_EXPR_FULL, &mx);
    int32_t *bp = malloc((size_t)(mx + 1) * 4), *kp = malloc((size_t)(mx + 1) * 4);
    int32_t *ij = malloc((size_t)(mx + 1) * 4), *pq = malloc((size_t)(mx + 1) * 4);
    double *cf = malloc((size_t)(mx + 1) * 8);
    if (hexpr_block_eval_full(drt, bra, nb, ket, nk, HEXPR_EXPR_FULL, &n, bp, kp, ij, pq, cf)
        != HEXPR_OK) FAIL("%s: block_eval_full status", name);
    hexpr_block_count_full(drt, bra, nb, ket, nk, HEXPR_EXPR_FULL, &cnt);
    if (cnt != n) FAIL("%s: block_count_full %d != block_eval_full %d", name, cnt, n);
    int k = 0;
    for (int i = 0; i < nb; i++)
        for (int j = 0; j < nk; j++) {
            tset_t P = get(drt, bra[i], ket[j], HEXPR_EXPR_FULL, 1);
            tset_t Q = {malloc((size_t)(P.n + 1) * sizeof(term_t)), 0};
            while (k < n && bp[k] == i && kp[k] == j) {
                Q.t[Q.n++] = (term_t){ij[k], pq[k], cf[k]};
                k++;
            }
            qsort(Q.t, (size_t)Q.n, sizeof *Q.t, term_cmp);
            if (!same(&P, &Q)) FAIL("%s: block (%d,%d) != pair_eval_full", name, i, j);
            free(P.t); free(Q.t);
        }
    if (k != n) FAIL("%s: block output not grouped by (bra_pos, ket_pos)", name);
    free(bp); free(kp); free(ij); free(pq); free(cf);
}

static unsigned rng = 12345u;
static int rnd(int n) { rng = rng * 1103515245u + 12345u; return (int)((rng >> 8) % (unsigned)n); }

static void run_drt(hexpr_drt_t *drt, const char *name, int max_all)
{
    if (!drt) { FAIL("%s: DRT build failed: %s", name, hexpr_last_error()); return; }
    int n = hexpr_drt_ncsfs(drt);
    long npairs = 0;
    if (n <= max_all) {
        for (int I = 0; I < n; I++)
            for (int J = I; J < n; J++) { check_pair(drt, name, I, J); npairs++; }
    } else {
        for (int t = 0; t < 3000; t++) { check_pair(drt, name, rnd(n), rnd(n)); npairs++; }
    }
    /* a disjoint off-diagonal block and a square one */
    int nb = n < 6 ? n / 2 : 3, nk = n < 6 ? n - nb : 4;
    int bra[3], ket[4], sq[4];
    for (int i = 0; i < nb; i++) bra[i] = i;
    for (int j = 0; j < nk; j++) ket[j] = n - 1 - j;
    for (int j = 0; j < 4; j++) sq[j] = (j * 7) % n;
    check_block(drt, name, bra, nb, ket, nk);
    check_block(drt, name, sq, n < 4 ? n : 4, sq, n < 4 ? n : 4);
    printf("  %-14s %4d CSFs, %6ld pairs checked\n", name, n, npairs);
}

int main(void)
{
    if (hexpr_init() != HEXPR_OK) { printf("FAIL: hexpr_init\n"); return 1; }
    printf("=== complete-decomposition calls (v1.1) vs union of directions ===\n");
    struct { const char *name; hexpr_drt_kind_t k; int nve, s2, nc, nv, ne; } cs[] = {
        {"s_full_4o4e", HEXPR_DRT_SOCI, 4, 0, 0, 4, 0},
        {"t_full_4o4e", HEXPR_DRT_SOCI, 4, 2, 0, 4, 0},
        {"q_full_4o4e", HEXPR_DRT_SOCI, 4, 4, 0, 4, 0},
        {"s_full_6o6e", HEXPR_DRT_SOCI, 6, 0, 0, 6, 0},
        {"t_full_6o6e", HEXPR_DRT_SOCI, 6, 2, 0, 6, 0},
        {"s_soci_small", HEXPR_DRT_SOCI, 6, 0, 1, 4, 1},
        {"t_soci_small", HEXPR_DRT_SOCI, 6, 2, 1, 4, 1},
        {"s_foci_small", HEXPR_DRT_FOCI, 6, 0, 1, 4, 1},
        {"t_foci_small", HEXPR_DRT_FOCI, 6, 2, 1, 4, 1},
        {"s_foci_ext2", HEXPR_DRT_FOCI, 6, 0, 1, 3, 3},
    };
    for (size_t c = 0; c < sizeof cs / sizeof cs[0]; c++) {
        hexpr_drt_t *d = hexpr_drt_create(cs[c].k, cs[c].nve, cs[c].s2, cs[c].nc, cs[c].nv, cs[c].ne);
        run_drt(d, cs[c].name, 200);
        if (d) hexpr_drt_destroy(d);
    }
    /* CSFSET sub-DRTs: scattered rows of a FOCI space, in shuffled order */
    hexpr_drt_t *base = hexpr_drt_create(HEXPR_DRT_FOCI, 6, 0, 1, 4, 1);
    if (base) {
        int n = hexpr_drt_ncsfs(base), norb = hexpr_drt_norb(base);
        const unsigned char *rows = hexpr_drt_csfs(base);
        for (int rep = 0; rep < 3; rep++) {
            int m = 25 + 10 * rep;
            unsigned char *buf = malloc((size_t)m * (size_t)norb);
            for (int i = 0; i < m; i++)
                memcpy(buf + (size_t)i * norb, rows + (size_t)rnd(n) * norb, (size_t)norb);
            hexpr_drt_t *sub = hexpr_drt_from_csfset(norb, buf, m);
            char name[32];
            snprintf(name, sizeof name, "csfset_%d", m);
            run_drt(sub, name, 200);
            if (sub) hexpr_drt_destroy(sub);
            free(buf);
        }
        hexpr_drt_destroy(base);
    }
    printf("  (pairs where the two v1.0 directions differed: %ld)\n", n_strict);
    if (failures) { printf("FAILED: %d check(s)\n", failures); return 1; }
    printf("ALL PASS\n");
    return 0;
}
