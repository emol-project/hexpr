/* test_pair_full_4o4e.c -- all 400 pairs vs pair.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "../include/hexpr.h"

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
/* Row type: (ij_key, pqrs, seq, coef)                                 */
/* seq records insertion order for first-wins dedup after sorting      */
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

/* Sort by (ij_key, pqrs, seq), then dedup keeping smallest seq (first-wins).
 * Returns new size. */
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

/* ------------------------------------------------------------------ */
/* Reference parsing                                                   */
/* ------------------------------------------------------------------ */
static int parse_reference(const char *path, row_vec_t *out)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return -1; }

    char line[512];
    /* skip to Energy Expression */
    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "Energy Expression")) break;
    }
    /* skip column header */
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }

    int I, J, pqrs, seq = 0;
    double coef;
    while (fscanf(f, "%d %d %d %lf", &I, &J, &pqrs, &coef) == 4) {
        rv_push(out, triangle(I-1, J-1), pqrs, seq++, coef);
    }
    fclose(f);
    return 0;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */
int main(void)
{
    printf("=== pair_eval full test (s_full_4o4e, all 400 pairs) ===\n");

    const char *ref_path = "../../reference_pair/s_full_4o4e/pair.txt";

    /* 1. Parse reference */
    row_vec_t ref = {NULL, 0, 0};
    if (parse_reference(ref_path, &ref) != 0) {
        printf("FAIL: cannot read reference file %s\n", ref_path);
        return 1;
    }
    int ref_raw = ref.size;
    int ref_n = sort_dedup(&ref);
    printf("Reference: %d raw entries -> %d unique (ij_key,pqrs) keys\n",
           ref_raw, ref_n);

    /* 2. Build DRT and evaluate all 400 pairs */
    hexpr_drt_t *drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 4, 0);
    if (!drt) { printf("FAIL: DRT creation failed\n"); return 1; }
    int ncsfs = hexpr_drt_ncsfs(drt);
    if (ncsfs != 20) {
        printf("FAIL: expected ncsfs=20, got %d\n", ncsfs);
        hexpr_drt_destroy(drt);
        return 1;
    }

    row_vec_t got = {NULL, 0, 0};
    int gseq = 0;
    for (int bra = 0; bra < ncsfs; bra++) {
        for (int ket = 0; ket < ncsfs; ket++) {
            int cnt = 0;
            hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_FULL, &cnt,
                            NULL, NULL, NULL);
            if (cnt == 0) continue;
            int32_t *ij_arr   = malloc((size_t)cnt * sizeof(int32_t));
            int32_t *pqrs_arr = malloc((size_t)cnt * sizeof(int32_t));
            double  *coef_arr = malloc((size_t)cnt * sizeof(double));
            hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_FULL, &cnt,
                            ij_arr, pqrs_arr, coef_arr);
            for (int t = 0; t < cnt; t++)
                rv_push(&got, ij_arr[t], pqrs_arr[t], gseq++, coef_arr[t]);
            free(ij_arr); free(pqrs_arr); free(coef_arr);
        }
    }
    int got_raw = got.size;
    int got_n = sort_dedup(&got);
    printf("C output:  %d raw entries -> %d unique (ij_key,pqrs) keys\n",
           got_raw, got_n);

    hexpr_drt_destroy(drt);

    /* 3. Compare */
    int key_missing_ref = 0;  /* in got but not ref */
    int key_missing_got = 0;  /* in ref but not got */
    int coef_mismatch   = 0;
    double max_adiff    = 0.0;
    int first_fail_printed = 0;

    int ri = 0, gi = 0;
    while (ri < ref_n || gi < got_n) {
        int cmp;
        if (ri >= ref_n) cmp = 1;
        else if (gi >= got_n) cmp = -1;
        else {
            if (ref.data[ri].ij_key != got.data[gi].ij_key)
                cmp = ref.data[ri].ij_key - got.data[gi].ij_key;
            else
                cmp = ref.data[ri].pqrs - got.data[gi].pqrs;
        }

        if (cmp < 0) {
            /* key in ref, not in got */
            if (!first_fail_printed) {
                printf("  FAIL: key (ij=%d,pqrs=%d) in ref but not C\n",
                       ref.data[ri].ij_key, ref.data[ri].pqrs);
                first_fail_printed = 1;
            }
            key_missing_got++;
            ri++;
        } else if (cmp > 0) {
            /* key in got, not in ref */
            if (!first_fail_printed) {
                printf("  FAIL: key (ij=%d,pqrs=%d) in C but not ref\n",
                       got.data[gi].ij_key, got.data[gi].pqrs);
                first_fail_printed = 1;
            }
            key_missing_ref++;
            gi++;
        } else {
            /* key present in both: compare coef */
            double adiff = fabs(ref.data[ri].coef - got.data[gi].coef);
            if (adiff > max_adiff) max_adiff = adiff;
            if (adiff > 1e-10) {
                coef_mismatch++;
                if (!first_fail_printed) {
                    printf("  FAIL: (ij=%d,pqrs=%d) ref=%.16e got=%.16e diff=%.3e\n",
                           ref.data[ri].ij_key, ref.data[ri].pqrs,
                           ref.data[ri].coef, got.data[gi].coef, adiff);
                    first_fail_printed = 1;
                }
            }
            ri++; gi++;
        }
    }

    free(ref.data);
    free(got.data);

    int failures = key_missing_ref + key_missing_got + coef_mismatch;
    printf("\nResults:\n");
    printf("  Keys in ref not in C: %d\n", key_missing_got);
    printf("  Keys in C not in ref: %d\n", key_missing_ref);
    printf("  Coef mismatches (>1e-10): %d\n", coef_mismatch);
    printf("  Max abs coef diff: %.3e\n", max_adiff);

    if (failures == 0)
        printf("\nALL TESTS PASSED\n");
    else
        printf("\n%d FAILURE(S)\n", failures);

    return failures > 0 ? 1 : 0;
}
