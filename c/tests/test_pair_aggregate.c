/* test_pair_aggregate.c -- generic pair aggregation test
 *
 * Usage: test_pair_aggregate <case_name> [<case_name> ...]
 *
 * For each case: evaluates hexpr_pair_eval for all (bra,ket) pairs,
 * aggregates with Hermitian first-wins deduplication (spec D7),
 * and compares against reference/<case>.expr.txt. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "../include/hexpr.h"

/* ------------------------------------------------------------------ */
/* Case table (mirrors run_reference.c ALL_CASES, 4o4e/6o6e only)     */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *name;
    int nve, spin, ncore, nval, next;
    hexpr_drt_kind_t kind;
} case_def_t;

static const case_def_t CASES[] = {
    {"s_full_4o4e", 4, 0, 0, 4, 0, HEXPR_DRT_SOCI},
    {"t_full_4o4e", 4, 2, 0, 4, 0, HEXPR_DRT_SOCI},
    {"s_full_6o6e", 6, 0, 0, 6, 0, HEXPR_DRT_SOCI},
    {"t_full_6o6e", 6, 2, 0, 6, 0, HEXPR_DRT_SOCI},
    {"s_foci_small", 6, 0, 1, 4, 1, HEXPR_DRT_FOCI},
    {"t_foci_small", 6, 2, 1, 4, 1, HEXPR_DRT_FOCI},
};
static const int NCASES = (int)(sizeof(CASES) / sizeof(CASES[0]));

/* ------------------------------------------------------------------ */
/* Row type and vector                                                 */
/* ------------------------------------------------------------------ */

typedef struct { int ij_key, pqrs, seq; double coef; } row_t;
typedef struct { row_t *data; int size, cap; } row_vec_t;

static void rv_push(row_vec_t *v, int ij_key, int pqrs, int seq, double coef)
{
    if (v->size == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 4096;
        v->data = realloc(v->data, (size_t)v->cap * sizeof(row_t));
    }
    v->data[v->size++] = (row_t){ij_key, pqrs, seq, coef};
}

static void rv_free(row_vec_t *v)
{
    free(v->data); v->data = NULL; v->size = v->cap = 0;
}

static int row_cmp(const void *a, const void *b)
{
    const row_t *ra = (const row_t *)a, *rb = (const row_t *)b;
    if (ra->ij_key != rb->ij_key) return ra->ij_key - rb->ij_key;
    if (ra->pqrs   != rb->pqrs)   return ra->pqrs   - rb->pqrs;
    return ra->seq - rb->seq;
}

/* Sort by (ij_key, pqrs, seq), dedup keeping first (smallest seq). Returns new size. */
static int sort_dedup(row_vec_t *v)
{
    if (v->size == 0) return 0;
    qsort(v->data, (size_t)v->size, sizeof(row_t), row_cmp);
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
/* Triangle index: triangle(a, b) = max*(max+1)/2 + min              */
/* ------------------------------------------------------------------ */

static int triangle(int a, int b)
{
    int mx = a > b ? a : b, mn = a < b ? a : b;
    return mx * (mx + 1) / 2 + mn;
}

/* ------------------------------------------------------------------ */
/* Reference parsing: ../../reference/<case>.expr.txt                  */
/* Format after "Energy Expression" + column header:                   */
/*   IJ  pqrs  coef   (IJ = triangle(bra,ket)+1, 1-based)             */
/* ------------------------------------------------------------------ */

static int parse_reference(const char *path, row_vec_t *out)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return -1; }

    char line[512];
    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, "Energy Expression")) break;
    }
    /* skip column header line */
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }

    int IJ, pqrs, seq = 0;
    double coef;
    while (fscanf(f, "%d %d %lf", &IJ, &pqrs, &coef) == 3) {
        rv_push(out, IJ - 1, pqrs, seq++, coef);
    }
    fclose(f);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Per-case runner                                                     */
/* ------------------------------------------------------------------ */

static int run_case(const case_def_t *c)
{
    printf("[%s] ... ", c->name);
    fflush(stdout);

    /* Build DRT */
    hexpr_drt_t *drt = hexpr_drt_create(c->kind, c->nve, c->spin,
                                        c->ncore, c->nval, c->next);
    if (!drt) { printf("FAIL: DRT creation failed\n"); return 1; }
    int ncsfs = hexpr_drt_ncsfs(drt);

    /* Evaluate all (bra, ket) pairs and aggregate */
    row_vec_t got = {NULL, 0, 0};
    int gseq = 0;
    int buf_cap = 65;
    int32_t *pqrs_buf = malloc((size_t)buf_cap * sizeof(int32_t));
    double  *coef_buf = malloc((size_t)buf_cap * sizeof(double));
    if (!pqrs_buf || !coef_buf) {
        printf("FAIL: OOM\n");
        hexpr_drt_destroy(drt);
        free(pqrs_buf); free(coef_buf);
        return 1;
    }

    for (int bra = 0; bra < ncsfs; bra++) {
        for (int ket = 0; ket < ncsfs; ket++) {
            int needed;
            hexpr_pair_count(drt, bra, ket, HEXPR_EXPR_FULL, &needed);
            if (needed == 0) continue;
            if (needed > buf_cap) {
                buf_cap = needed * 2;
                pqrs_buf = realloc(pqrs_buf, (size_t)buf_cap * sizeof(int32_t));
                coef_buf = realloc(coef_buf, (size_t)buf_cap * sizeof(double));
            }
            int cnt;
            hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_FULL,
                            &cnt, NULL, pqrs_buf, coef_buf);
            int ij_key = triangle(bra, ket);
            for (int t = 0; t < cnt; t++)
                rv_push(&got, ij_key, (int)pqrs_buf[t], gseq++, coef_buf[t]);
        }
    }

    free(pqrs_buf);
    free(coef_buf);
    hexpr_drt_destroy(drt);

    int got_n = sort_dedup(&got);

    /* Parse reference */
    char ref_path[512];
    snprintf(ref_path, sizeof(ref_path), "../../reference/%s.expr.txt", c->name);
    row_vec_t ref = {NULL, 0, 0};
    if (parse_reference(ref_path, &ref) != 0) {
        printf("FAIL: cannot read reference\n");
        rv_free(&got);
        return 1;
    }
    int ref_n = sort_dedup(&ref);

    /* Compare */
    int key_missing_ref = 0;
    int key_missing_got = 0;
    int coef_mismatch   = 0;
    double max_adiff    = 0.0;
    int n_reported      = 0;

    int ri = 0, gi = 0;
    while (ri < ref_n || gi < got_n) {
        int cmp;
        if      (ri >= ref_n) cmp =  1;
        else if (gi >= got_n) cmp = -1;
        else if (ref.data[ri].ij_key != got.data[gi].ij_key)
            cmp = ref.data[ri].ij_key - got.data[gi].ij_key;
        else
            cmp = ref.data[ri].pqrs - got.data[gi].pqrs;

        if (cmp < 0) {
            if (n_reported++ < 3)
                fprintf(stderr, "  FAIL [%s]: (ij=%d,pqrs=%d) in ref but not C\n",
                        c->name, ref.data[ri].ij_key, ref.data[ri].pqrs);
            key_missing_got++;
            ri++;
        } else if (cmp > 0) {
            if (n_reported++ < 3)
                fprintf(stderr, "  FAIL [%s]: (ij=%d,pqrs=%d) in C but not ref\n",
                        c->name, got.data[gi].ij_key, got.data[gi].pqrs);
            key_missing_ref++;
            gi++;
        } else {
            double adiff = fabs(ref.data[ri].coef - got.data[gi].coef);
            if (adiff > max_adiff) max_adiff = adiff;
            if (adiff > 1e-10) {
                if (n_reported++ < 3)
                    fprintf(stderr,
                            "  FAIL [%s]: (ij=%d,pqrs=%d) ref=%.16e got=%.16e diff=%.3e\n",
                            c->name, ref.data[ri].ij_key, ref.data[ri].pqrs,
                            ref.data[ri].coef, got.data[gi].coef, adiff);
                coef_mismatch++;
            }
            ri++; gi++;
        }
    }

    rv_free(&ref);
    rv_free(&got);

    int failures = key_missing_ref + key_missing_got + coef_mismatch;
    if (failures == 0) {
        printf("PASS (ref=%d C=%d max_diff=%.3e)\n", ref_n, got_n, max_adiff);
    } else {
        printf("FAIL (missing_in_C=%d missing_in_ref=%d coef_mismatch=%d max_diff=%.3e)\n",
               key_missing_got, key_missing_ref, coef_mismatch, max_adiff);
    }
    return failures > 0 ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <case_name> [<case_name> ...]\n", argv[0]);
        fprintf(stderr, "Known cases:");
        for (int i = 0; i < NCASES; i++) fprintf(stderr, " %s", CASES[i].name);
        fprintf(stderr, "\n");
        return 1;
    }

    int total_fail = 0;
    for (int a = 1; a < argc; a++) {
        const case_def_t *c = NULL;
        for (int i = 0; i < NCASES; i++)
            if (strcmp(CASES[i].name, argv[a]) == 0) { c = &CASES[i]; break; }
        if (!c) {
            fprintf(stderr, "Unknown case: %s\n", argv[a]);
            total_fail++;
            continue;
        }
        if (run_case(c) != 0) total_fail++;
    }

    if (total_fail == 0)
        printf("ALL PASS\n");
    else
        printf("%d FAIL(S)\n", total_fail);

    return total_fail > 0 ? 1 : 0;
}
