/* test_block_subspace_4o4e.c -- block_eval / subspace_eval */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "../include/hexpr.h"

static int failures = 0;
static int checks   = 0;

static void check(const char *label, int ok)
{
    checks++;
    printf("  %-62s %s\n", label, ok ? "PASS" : "FAIL");
    if (!ok) failures++;
}

static void check_status(const char *label, hexpr_status_t got,
                          hexpr_status_t expected)
{
    checks++;
    int ok = (got == expected);
    printf("  %-62s %s  (got %d, expected %d)\n",
           label, ok ? "PASS" : "FAIL", (int)got, (int)expected);
    if (!ok) failures++;
}

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
/* ------------------------------------------------------------------ */
static int imax2(int a, int b) { return a > b ? a : b; }
static int imin2(int a, int b) { return a < b ? a : b; }
static int triangle(int a, int b)
{
    int mx = imax2(a,b), mn = imin2(a,b);
    return mx*(mx+1)/2 + mn;
}

/* ------------------------------------------------------------------ */
/* First-wins dedup helpers (same as test_pair_full_4o4e)              */
/* ------------------------------------------------------------------ */
typedef struct { int ij_key, pqrs, seq; double coef; } row_t;
typedef struct { row_t *data; int size, cap; } row_vec_t;

static void rv_push(row_vec_t *v, int ij_key, int pqrs, int seq, double coef)
{
    if (v->size == v->cap) {
        v->cap = v->cap ? v->cap*2 : 4096;
        v->data = realloc(v->data, v->cap * sizeof(row_t));
    }
    v->data[v->size++] = (row_t){ij_key, pqrs, seq, coef};
}

static int row_cmp(const void *a, const void *b)
{
    const row_t *ra = (const row_t *)a, *rb = (const row_t *)b;
    if (ra->ij_key != rb->ij_key) return ra->ij_key - rb->ij_key;
    if (ra->pqrs   != rb->pqrs)   return ra->pqrs   - rb->pqrs;
    return ra->seq - rb->seq;
}

static int sort_dedup(row_vec_t *v)
{
    if (v->size == 0) return 0;
    qsort(v->data, v->size, sizeof(row_t), row_cmp);
    int out = 0;
    for (int i = 0; i < v->size; ) {
        v->data[out++] = v->data[i];
        int j = i + 1;
        while (j < v->size &&
               v->data[j].ij_key == v->data[i].ij_key &&
               v->data[j].pqrs   == v->data[i].pqrs) j++;
        i = j;
    }
    v->size = out;
    return out;
}

static int parse_reference(const char *path, row_vec_t *out)
{
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    char line[512];
    while (fgets(line, sizeof(line), f))
        if (strstr(line, "Energy Expression")) break;
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }
    int I, J, pqrs, seq = 0;
    double coef;
    while (fscanf(f, "%d %d %d %lf", &I, &J, &pqrs, &coef) == 4)
        rv_push(out, triangle(I-1, J-1), pqrs, seq++, coef);
    fclose(f);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Test 1: block_eval 1x1 == pair_eval for sampled pairs               */
/* ------------------------------------------------------------------ */
static void test1_block_1x1(hexpr_drt_t *drt)
{
    printf("\n--- Test 1: block_eval 1x1 == pair_eval (sampled pairs) ---\n");

    int pairs[][2] = {{0,0},{0,5},{7,3},{19,19}};
    int npairs = (int)(sizeof(pairs)/sizeof(pairs[0]));

    for (int p = 0; p < npairs; p++) {
        int bra = pairs[p][0], ket = pairs[p][1];

        /* pair_eval */
        int pcnt = 0;
        hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_FULL, &pcnt,
                        NULL, NULL, NULL);
        int32_t *pij   = malloc((size_t)(pcnt+1) * sizeof(int32_t));
        int32_t *ppqrs = malloc((size_t)(pcnt+1) * sizeof(int32_t));
        double  *pcoef = malloc((size_t)(pcnt+1) * sizeof(double));
        hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_FULL, &pcnt,
                        pij, ppqrs, pcoef);

        /* block_eval 1x1 */
        int bcnt = 0;
        int bra_arr[1] = {bra}, ket_arr[1] = {ket};
        hexpr_block_eval(drt, bra_arr, 1, ket_arr, 1, HEXPR_EXPR_FULL,
                         &bcnt, NULL, NULL, NULL, NULL, NULL);

        char label[128];
        snprintf(label, sizeof(label),
                 "block_count 1x1 (%d,%d) == pair_count", bra, ket);
        check(label, bcnt == pcnt);

        if (bcnt == pcnt && bcnt > 0) {
            int32_t *bij   = malloc((size_t)bcnt * sizeof(int32_t));
            int32_t *bpqrs = malloc((size_t)bcnt * sizeof(int32_t));
            double  *bcoef = malloc((size_t)bcnt * sizeof(double));
            int32_t *bbra  = malloc((size_t)bcnt * sizeof(int32_t));
            int32_t *bket  = malloc((size_t)bcnt * sizeof(int32_t));
            hexpr_block_eval(drt, bra_arr, 1, ket_arr, 1, HEXPR_EXPR_FULL,
                             &bcnt, bbra, bket, bij, bpqrs, bcoef);

            int all_ok = 1;
            for (int t = 0; t < bcnt; t++) {
                if (bij[t] != pij[t] || bpqrs[t] != ppqrs[t] ||
                    fabs(bcoef[t] - pcoef[t]) > 1e-14 ||
                    bbra[t] != 0 || bket[t] != 0) {
                    all_ok = 0; break;
                }
            }
            snprintf(label, sizeof(label),
                     "block_eval 1x1 (%d,%d) terms == pair_eval terms", bra, ket);
            check(label, all_ok);

            free(bij); free(bpqrs); free(bcoef); free(bbra); free(bket);
        }

        free(pij); free(ppqrs); free(pcoef);
    }
}

/* ------------------------------------------------------------------ */
/* Test 2: block_eval rectangular 2x3                                  */
/* ------------------------------------------------------------------ */
static void test2_block_2x3(hexpr_drt_t *drt)
{
    printf("\n--- Test 2: block_eval rectangular 2x3 ---\n");

    int bra_arr[2] = {0, 5};
    int ket_arr[3] = {1, 7, 12};
    int n_bra = 2, n_ket = 3;

    /* Expected total = sum of pair counts */
    int expected_total = 0;
    for (int ip = 0; ip < n_bra; ip++) {
        for (int jp = 0; jp < n_ket; jp++) {
            int cnt = 0;
            hexpr_pair_eval(drt, bra_arr[ip], ket_arr[jp], HEXPR_EXPR_FULL,
                            &cnt, NULL, NULL, NULL);
            expected_total += cnt;
        }
    }

    int bcnt = 0;
    hexpr_block_count(drt, bra_arr, n_bra, ket_arr, n_ket,
                      HEXPR_EXPR_FULL, &bcnt);
    check("block_count 2x3 == sum of pair counts", bcnt == expected_total);

    if (bcnt != expected_total) return;

    int32_t *bbra  = malloc((size_t)bcnt * sizeof(int32_t));
    int32_t *bket  = malloc((size_t)bcnt * sizeof(int32_t));
    int32_t *bij   = malloc((size_t)bcnt * sizeof(int32_t));
    int32_t *bpqrs = malloc((size_t)bcnt * sizeof(int32_t));
    double  *bcoef = malloc((size_t)bcnt * sizeof(double));

    hexpr_block_eval(drt, bra_arr, n_bra, ket_arr, n_ket, HEXPR_EXPR_FULL,
                     &bcnt, bbra, bket, bij, bpqrs, bcoef);

    /* Verify term tagging: group terms by (bra_pos, ket_pos) and compare
     * to pair_eval for the corresponding pair. */
    int all_ok = 1;
    int t = 0;
    for (int ip = 0; ip < n_bra && all_ok; ip++) {
        for (int jp = 0; jp < n_ket && all_ok; jp++) {
            int pcnt = 0;
            hexpr_pair_eval(drt, bra_arr[ip], ket_arr[jp], HEXPR_EXPR_FULL,
                            &pcnt, NULL, NULL, NULL);
            int32_t *pij   = malloc((size_t)(pcnt+1) * sizeof(int32_t));
            int32_t *ppqrs = malloc((size_t)(pcnt+1) * sizeof(int32_t));
            double  *pcoef = malloc((size_t)(pcnt+1) * sizeof(double));
            hexpr_pair_eval(drt, bra_arr[ip], ket_arr[jp], HEXPR_EXPR_FULL,
                            &pcnt, pij, ppqrs, pcoef);

            for (int k = 0; k < pcnt; k++, t++) {
                if (t >= bcnt || bbra[t] != ip || bket[t] != jp ||
                    bij[t] != pij[k] || bpqrs[t] != ppqrs[k] ||
                    fabs(bcoef[t] - pcoef[k]) > 1e-14) {
                    all_ok = 0; break;
                }
            }
            free(pij); free(ppqrs); free(pcoef);
        }
    }
    check("block_eval 2x3 term tagging matches pair_eval", all_ok);

    free(bbra); free(bket); free(bij); free(bpqrs); free(bcoef);
}

/* ------------------------------------------------------------------ */
/* Test 3: subspace_eval == block_eval(indices, indices)               */
/* ------------------------------------------------------------------ */
static void test3_subspace_eq_block(hexpr_drt_t *drt)
{
    printf("\n--- Test 3: subspace_eval == block_eval(indices, indices) ---\n");

    int indices[3] = {0, 1, 2};
    int n = 3;

    /* block_eval with same arrays */
    int bcnt = 0;
    hexpr_block_count(drt, indices, n, indices, n, HEXPR_EXPR_FULL, &bcnt);

    int scnt = 0;
    hexpr_subspace_count(drt, indices, n, HEXPR_EXPR_FULL, &scnt);
    check("subspace_count([0,1,2]) == block_count([0,1,2],[0,1,2])",
          scnt == bcnt);

    if (scnt != bcnt) return;

    int32_t *bbra  = malloc((size_t)bcnt * sizeof(int32_t));
    int32_t *bket  = malloc((size_t)bcnt * sizeof(int32_t));
    int32_t *bij   = malloc((size_t)bcnt * sizeof(int32_t));
    int32_t *bpqrs = malloc((size_t)bcnt * sizeof(int32_t));
    double  *bcoef = malloc((size_t)bcnt * sizeof(double));
    hexpr_block_eval(drt, indices, n, indices, n, HEXPR_EXPR_FULL,
                     &bcnt, bbra, bket, bij, bpqrs, bcoef);

    int32_t *sbra  = malloc((size_t)scnt * sizeof(int32_t));
    int32_t *sket  = malloc((size_t)scnt * sizeof(int32_t));
    int32_t *sij   = malloc((size_t)scnt * sizeof(int32_t));
    int32_t *spqrs = malloc((size_t)scnt * sizeof(int32_t));
    double  *scoef = malloc((size_t)scnt * sizeof(double));
    hexpr_subspace_eval(drt, indices, n, HEXPR_EXPR_FULL,
                        &scnt, sbra, sket, sij, spqrs, scoef);

    int ok = (memcmp(bbra,  sbra,  bcnt*sizeof(int32_t)) == 0 &&
              memcmp(bket,  sket,  bcnt*sizeof(int32_t)) == 0 &&
              memcmp(bij,   sij,   bcnt*sizeof(int32_t)) == 0 &&
              memcmp(bpqrs, spqrs, bcnt*sizeof(int32_t)) == 0 &&
              memcmp(bcoef, scoef, bcnt*sizeof(double))  == 0);
    check("subspace_eval([0,1,2]) byte-identical to block_eval", ok);

    free(bbra); free(bket); free(bij); free(bpqrs); free(bcoef);
    free(sbra); free(sket); free(sij); free(spqrs); free(scoef);
}

/* ------------------------------------------------------------------ */
/* Test 4: subspace_eval full (20 CSFs) == reference pair.txt          */
/* ------------------------------------------------------------------ */
static void test4_subspace_full_vs_reference(hexpr_drt_t *drt)
{
    printf("\n--- Test 4: subspace_eval full (20 CSFs) == pair.txt reference ---\n");

    const char *ref_path = "../../reference_pair/s_full_4o4e/pair.txt";
    row_vec_t ref = {NULL, 0, 0};
    if (parse_reference(ref_path, &ref) != 0) {
        check("parse reference pair.txt", 0);
        return;
    }
    int ref_n = sort_dedup(&ref);

    /* Build index array [0..19] */
    int ncsfs = hexpr_drt_ncsfs(drt);
    int *indices = malloc((size_t)ncsfs * sizeof(int));
    for (int i = 0; i < ncsfs; i++) indices[i] = i;

    int scnt = 0;
    hexpr_subspace_eval(drt, indices, ncsfs, HEXPR_EXPR_FULL,
                        &scnt, NULL, NULL, NULL, NULL, NULL);

    int32_t *sbra  = malloc((size_t)scnt * sizeof(int32_t));
    int32_t *sket  = malloc((size_t)scnt * sizeof(int32_t));
    int32_t *sij   = malloc((size_t)scnt * sizeof(int32_t));
    int32_t *spqrs = malloc((size_t)scnt * sizeof(int32_t));
    double  *scoef = malloc((size_t)scnt * sizeof(double));
    hexpr_subspace_eval(drt, indices, ncsfs, HEXPR_EXPR_FULL,
                        &scnt, sbra, sket, sij, spqrs, scoef);

    /* Aggregate with first-wins */
    row_vec_t got = {NULL, 0, 0};
    for (int t = 0; t < scnt; t++)
        rv_push(&got, sij[t], spqrs[t], t, scoef[t]);
    int got_n = sort_dedup(&got);

    free(sbra); free(sket); free(sij); free(spqrs); free(scoef);
    free(indices);

    check("unique key count matches reference", got_n == ref_n);

    int mismatches = 0;
    double max_diff = 0.0;
    int ri = 0, gi = 0;
    while (ri < ref_n && gi < got_n) {
        int cmp;
        if (ref.data[ri].ij_key != got.data[gi].ij_key)
            cmp = ref.data[ri].ij_key - got.data[gi].ij_key;
        else
            cmp = ref.data[ri].pqrs - got.data[gi].pqrs;

        if (cmp < 0) { mismatches++; ri++; }
        else if (cmp > 0) { mismatches++; gi++; }
        else {
            double d = fabs(ref.data[ri].coef - got.data[gi].coef);
            if (d > max_diff) max_diff = d;
            if (d > 1e-10) mismatches++;
            ri++; gi++;
        }
    }
    mismatches += (ref_n - ri) + (got_n - gi);

    printf("  Max abs coef diff vs reference: %.3e\n", max_diff);
    check("zero mismatches vs pair.txt (Hermitian first-wins)", mismatches == 0);

    free(ref.data);
    free(got.data);
}

/* ------------------------------------------------------------------ */
/* Test 5: HEXPR_EXPR_ONE filter on block and subspace                 */
/* ------------------------------------------------------------------ */
static void test5_one_filter(hexpr_drt_t *drt)
{
    printf("\n--- Test 5: HEXPR_EXPR_ONE filter ---\n");

    int norb  = hexpr_drt_norb(drt);
    int n_oel = norb * (norb + 1) / 2;

    int bra_arr[2] = {0, 3};
    int ket_arr[2] = {0, 3};

    int cnt_full = 0, cnt_one = 0;
    hexpr_block_count(drt, bra_arr, 2, ket_arr, 2,
                      HEXPR_EXPR_FULL, &cnt_full);
    hexpr_block_count(drt, bra_arr, 2, ket_arr, 2,
                      HEXPR_EXPR_ONE, &cnt_one);
    check("block: ONE count <= FULL count", cnt_one <= cnt_full);

    if (cnt_one > 0) {
        int32_t *pqrs_arr = malloc((size_t)cnt_one * sizeof(int32_t));
        hexpr_block_eval(drt, bra_arr, 2, ket_arr, 2, HEXPR_EXPR_ONE,
                         &cnt_one, NULL, NULL, NULL, pqrs_arr, NULL);
        int all_1e = 1;
        for (int t = 0; t < cnt_one; t++)
            if (pqrs_arr[t] >= n_oel) { all_1e = 0; break; }
        check("block ONE: all pqrs < n_oel", all_1e);
        free(pqrs_arr);
    } else {
        check("block ONE: (no terms to check)", 1);
    }

    int indices[3] = {0, 5, 10};
    int scnt_full = 0, scnt_one = 0;
    hexpr_subspace_count(drt, indices, 3, HEXPR_EXPR_FULL, &scnt_full);
    hexpr_subspace_count(drt, indices, 3, HEXPR_EXPR_ONE,  &scnt_one);
    check("subspace: ONE count <= FULL count", scnt_one <= scnt_full);

    if (scnt_one > 0) {
        int32_t *pqrs_arr = malloc((size_t)scnt_one * sizeof(int32_t));
        hexpr_subspace_eval(drt, indices, 3, HEXPR_EXPR_ONE,
                            &scnt_one, NULL, NULL, NULL, pqrs_arr, NULL);
        int all_1e = 1;
        for (int t = 0; t < scnt_one; t++)
            if (pqrs_arr[t] >= n_oel) { all_1e = 0; break; }
        check("subspace ONE: all pqrs < n_oel", all_1e);
        free(pqrs_arr);
    } else {
        check("subspace ONE: (no terms to check)", 1);
    }
}

/* ------------------------------------------------------------------ */
/* Test 6: invalid inputs                                              */
/* ------------------------------------------------------------------ */
static void test6_invalid(hexpr_drt_t *drt)
{
    printf("\n--- Test 6: invalid inputs ---\n");

    int ncsfs = hexpr_drt_ncsfs(drt);
    int ok_arr[2] = {0, 1};
    int bad_neg[2] = {0, -1};
    int bad_hi[2]  = {0, ncsfs};
    int cnt;

    /* NULL drt */
    check_status("block_count NULL drt -> INVALID_ARG",
        hexpr_block_count(NULL, ok_arr, 1, ok_arr, 1, HEXPR_EXPR_FULL, &cnt),
        HEXPR_ERR_INVALID_ARG);

    /* NULL out_count */
    check_status("block_count NULL out_count -> INVALID_ARG",
        hexpr_block_count(drt, ok_arr, 1, ok_arr, 1, HEXPR_EXPR_FULL, NULL),
        HEXPR_ERR_INVALID_ARG);

    /* n_bra < 0 */
    check_status("block_count n_bra=-1 -> INVALID_ARG",
        hexpr_block_count(drt, ok_arr, -1, ok_arr, 1, HEXPR_EXPR_FULL, &cnt),
        HEXPR_ERR_INVALID_ARG);

    /* bra index out of range */
    check_status("block_count bra_indices[1]=-1 -> RANGE",
        hexpr_block_count(drt, bad_neg, 2, ok_arr, 1, HEXPR_EXPR_FULL, &cnt),
        HEXPR_ERR_RANGE);
    check_status("block_count bra_indices[1]=ncsfs -> RANGE",
        hexpr_block_count(drt, bad_hi, 2, ok_arr, 1, HEXPR_EXPR_FULL, &cnt),
        HEXPR_ERR_RANGE);

    /* ket index out of range */
    check_status("block_count ket_indices[1]=-1 -> RANGE",
        hexpr_block_count(drt, ok_arr, 1, bad_neg, 2, HEXPR_EXPR_FULL, &cnt),
        HEXPR_ERR_RANGE);

    /* subspace NULL drt */
    check_status("subspace_count NULL drt -> INVALID_ARG",
        hexpr_subspace_count(NULL, ok_arr, 1, HEXPR_EXPR_FULL, &cnt),
        HEXPR_ERR_INVALID_ARG);

    /* subspace index out of range */
    check_status("subspace_count indices[1]=-1 -> RANGE",
        hexpr_subspace_count(drt, bad_neg, 2, HEXPR_EXPR_FULL, &cnt),
        HEXPR_ERR_RANGE);
}

/* ------------------------------------------------------------------ */
/* Test 7: zero-size inputs                                            */
/* ------------------------------------------------------------------ */
static void test7_zero_size(hexpr_drt_t *drt)
{
    printf("\n--- Test 7: zero-size inputs ---\n");

    int count;
    int dummy_idx[1] = {0};

    check("block_count n_bra=0 returns 0",
        hexpr_block_count(drt, NULL, 0, dummy_idx, 1,
                           HEXPR_EXPR_FULL, &count) == HEXPR_OK
        && count == 0);
    check("block_count n_ket=0 returns 0",
        hexpr_block_count(drt, dummy_idx, 1, NULL, 0,
                           HEXPR_EXPR_FULL, &count) == HEXPR_OK
        && count == 0);
    check("block_count n_bra=0 AND n_ket=0 returns 0",
        hexpr_block_count(drt, NULL, 0, NULL, 0,
                           HEXPR_EXPR_FULL, &count) == HEXPR_OK
        && count == 0);
    check("block_eval n_bra=0 returns 0",
        hexpr_block_eval(drt, NULL, 0, dummy_idx, 1,
                          HEXPR_EXPR_FULL, &count,
                          NULL, NULL, NULL, NULL, NULL) == HEXPR_OK
        && count == 0);
    check("subspace_count n=0 returns 0",
        hexpr_subspace_count(drt, NULL, 0,
                              HEXPR_EXPR_FULL, &count) == HEXPR_OK
        && count == 0);
    check("subspace_eval n=0 returns 0",
        hexpr_subspace_eval(drt, NULL, 0,
                             HEXPR_EXPR_FULL, &count,
                             NULL, NULL, NULL, NULL, NULL) == HEXPR_OK
        && count == 0);
}

/* ------------------------------------------------------------------ */
/* Test 8: duplicated indices                                           */
/* ------------------------------------------------------------------ */
static void test8_duplicated_indices(hexpr_drt_t *drt)
{
    printf("\n--- Test 8: duplicated indices ---\n");

    /* block_eval with bra=[0,0], ket=[1] should give 2 pair evaluations
     * of (0,1). The output should contain 2x the terms of pair_eval(0,1). */
    int n_pair, n_block;
    hexpr_pair_count(drt, 0, 1, HEXPR_EXPR_FULL, &n_pair);
    int bra_dup[2] = {0, 0};
    int ket[1] = {1};
    hexpr_block_count(drt, bra_dup, 2, ket, 1, HEXPR_EXPR_FULL, &n_block);
    check("block_count([0,0],[1]) == 2 * pair_count(0,1)",
          n_block == 2 * n_pair);

    /* subspace_count with duplicated indices */
    int n_pair_00, n_pair_01, n_pair_10, n_pair_11;
    hexpr_pair_count(drt, 0, 0, HEXPR_EXPR_FULL, &n_pair_00);
    hexpr_pair_count(drt, 0, 1, HEXPR_EXPR_FULL, &n_pair_01);
    hexpr_pair_count(drt, 1, 0, HEXPR_EXPR_FULL, &n_pair_10);
    hexpr_pair_count(drt, 1, 1, HEXPR_EXPR_FULL, &n_pair_11);
    int indices_dup[3] = {0, 0, 1};
    int n_sub;
    hexpr_subspace_count(drt, indices_dup, 3, HEXPR_EXPR_FULL, &n_sub);
    /* subspace([0,0,1]) = 3x3 = 9 pairs:
     *   (0,0) (0,0) (0,1)
     *   (0,0) (0,0) (0,1)
     *   (1,0) (1,0) (1,1)
     * = 4*pair(0,0) + 2*pair(0,1) + 2*pair(1,0) + 1*pair(1,1) */
    int expected = 4 * n_pair_00 + 2 * n_pair_01 + 2 * n_pair_10 + n_pair_11;
    check("subspace_count([0,0,1]) handles duplicates correctly",
          n_sub == expected);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */
int main(void)
{
    printf("=== block_eval / subspace_eval test (s_full_4o4e) ===\n");

    hexpr_drt_t *drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 4, 0);
    if (!drt) { printf("FAIL: DRT creation failed\n"); return 1; }
    check("DRT create SOCI 4o4e", 1);
    check("ncsfs == 20", hexpr_drt_ncsfs(drt) == 20);

    test1_block_1x1(drt);
    test2_block_2x3(drt);
    test3_subspace_eq_block(drt);
    test4_subspace_full_vs_reference(drt);
    test5_one_filter(drt);
    test6_invalid(drt);
    test7_zero_size(drt);
    test8_duplicated_indices(drt);

    hexpr_drt_destroy(drt);

    printf("\n=== block/subspace: %d/%d checks passed ===\n",
           checks - failures, checks);
    if (failures == 0)
        printf("ALL TESTS PASSED\n");
    else
        printf("%d TEST(S) FAILED\n", failures);
    return failures > 0 ? 1 : 0;
}
