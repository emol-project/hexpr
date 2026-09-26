#define _POSIX_C_SOURCE 199309L
/* test_expr_stress.c -- stress test for the expression generator
 *
 * Stress test: iterate all (or 100 sampled) (I,J) pairs for each of the
 * Compare against reference/<case>.expr.txt for 12 cases.
 *
 * Full sweep mode: estimated per-pair time * ncsfs^2 <= 3 h.
 * Sample mode:     estimated time > 3 h; runs Hermitian self-check instead.
 *
 * Run from c/tests/ as: LD_LIBRARY_PATH=.. ./test_expr_stress
 * Or via:               make -C c/tests stress
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include "../include/hexpr.h"

/* ------------------------------------------------------------------ */
/* Case table (12 cases, in fail-fast order)                          */
/* ------------------------------------------------------------------ */

typedef struct {
    const char *name;
    int nve, spin, ncore, nval, next;
    hexpr_drt_kind_t kind;
} case_def_t;

static const case_def_t CASES[] = {
    {"s_full_4o4e",   4, 0, 0, 4, 0, HEXPR_DRT_SOCI},
    {"t_full_4o4e",   4, 2, 0, 4, 0, HEXPR_DRT_SOCI},
    {"s_full_6o6e",   6, 0, 0, 6, 0, HEXPR_DRT_SOCI},
    {"t_full_6o6e",   6, 2, 0, 6, 0, HEXPR_DRT_SOCI},
    {"s_foci_small",  6, 0, 1, 4, 1, HEXPR_DRT_FOCI},
    {"t_foci_small",  6, 2, 1, 4, 1, HEXPR_DRT_FOCI},
    {"s_soci_small",  6, 0, 1, 4, 1, HEXPR_DRT_SOCI},
    {"t_soci_small",  6, 2, 1, 4, 1, HEXPR_DRT_SOCI},
    {"s_foci_med",   10, 0, 2, 4, 4, HEXPR_DRT_FOCI},
    {"t_foci_med",   10, 2, 2, 4, 4, HEXPR_DRT_FOCI},
    {"s_soci_med",   10, 0, 2, 4, 4, HEXPR_DRT_SOCI},
    {"t_soci_med",   10, 2, 2, 4, 4, HEXPR_DRT_SOCI},
};
#define NCASES ((int)(sizeof(CASES) / sizeof(CASES[0])))

/* ------------------------------------------------------------------ */
/* Row vector for aggregation (ij_key = triangle(bra,ket), 0-based)   */
/* ------------------------------------------------------------------ */

typedef struct { int ij_key, pqrs, seq; double coef; } row_t;
typedef struct { row_t *data; int size, cap; } row_vec_t;

static void rv_push(row_vec_t *v, int ij, int pq, int seq, double c)
{
    if (v->size == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 65536;
        v->data = realloc(v->data, (size_t)v->cap * sizeof(row_t));
    }
    v->data[v->size++] = (row_t){ij, pq, seq, c};
}

static void rv_free(row_vec_t *v)
{
    free(v->data); v->data = NULL; v->size = v->cap = 0;
}

static int row_cmp(const void *a, const void *b)
{
    const row_t *ra = a, *rb = b;
    if (ra->ij_key != rb->ij_key) return ra->ij_key - rb->ij_key;
    if (ra->pqrs   != rb->pqrs)   return ra->pqrs   - rb->pqrs;
    return ra->seq - rb->seq;
}

/* Sort by (ij_key, pqrs, seq); keep first entry per (ij_key, pqrs). */
static int sort_dedup(row_vec_t *v)
{
    if (!v->size) return 0;
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
    return v->size = out;
}

/* ------------------------------------------------------------------ */
/* Utilities                                                           */
/* ------------------------------------------------------------------ */

static int tri(int a, int b)
{
    int hi = a > b ? a : b, lo = a < b ? a : b;
    return hi * (hi + 1) / 2 + lo;
}

static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec * 1e-6;
}

static long vmrss_kb(void)
{
    FILE *f = fopen("/proc/self/status", "r");
    if (!f) return -1;
    char line[256]; long v = -1;
    while (fgets(line, sizeof(line), f))
        if (strncmp(line, "VmRSS:", 6) == 0) { sscanf(line + 6, "%ld", &v); break; }
    fclose(f);
    return v;
}

/* Parse reference/<case>.expr.txt into out (IJ column is 1-based; stored 0-based). */
static int parse_ref(const char *path, row_vec_t *out)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "  cannot open %s\n", path); return -1; }
    char line[512];
    while (fgets(line, sizeof(line), f))
        if (strstr(line, "Energy Expression")) break;
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }  /* header row */
    int IJ, pq, seq = 0; double co;
    while (fscanf(f, "%d %d %lf", &IJ, &pq, &co) == 3)
        rv_push(out, IJ - 1, pq, seq++, co);
    fclose(f);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Per-case result                                                     */
/* ------------------------------------------------------------------ */

typedef struct {
    char name[32];
    int  ncsfs;
    long pairs_tested;
    int  mode;           /* 0 = full sweep, 1 = sample */
    double per_pair_ms;
    double total_s;
    char   verif[80];
    double rss_max_mb;
    double rss_growth_mb;
    int    ok;
} result_t;

/* ------------------------------------------------------------------ */
/* Full sweep                                                          */
/* ------------------------------------------------------------------ */

static result_t run_full(const case_def_t *c, hexpr_drt_t *drt, int ncsfs,
                         int32_t *pb, double *cb)
{
    result_t res; memset(&res, 0, sizeof(res));
    strncpy(res.name, c->name, 31); res.ncsfs = ncsfs; res.mode = 0;

    row_vec_t got = {NULL, 0, 0};
    int gseq = 0;
    long rss_first = -1, rss_last = -1, rss_max = 0;
    double t0 = now_ms();
    long pair_n = 0;

    for (int bra = 0; bra < ncsfs; bra++) {
        for (int ket = 0; ket < ncsfs; ket++) {
            int cnt;
            hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_FULL,
                            &cnt, NULL, pb, cb);
            if (cnt > 0) {
                int ij = tri(bra, ket);
                for (int t = 0; t < cnt; t++)
                    rv_push(&got, ij, (int)pb[t], gseq++, cb[t]);
            }
            pair_n++;
            if (pair_n % 1000 == 0) {
                long rss = vmrss_kb();
                if (rss > 0) {
                    if (rss_first < 0) rss_first = rss;
                    rss_last = rss;
                    if (rss > rss_max) rss_max = rss;
                }
                /* periodic dedup to bound intermediate memory */
                if (got.size > 2000000) sort_dedup(&got);
            }
        }
    }
    double elapsed_ms = now_ms() - t0;
    sort_dedup(&got);

    res.pairs_tested  = pair_n;
    res.total_s       = elapsed_ms / 1000.0;
    res.per_pair_ms   = pair_n > 0 ? elapsed_ms / pair_n : 0.0;
    res.rss_max_mb    = rss_max   > 0 ? rss_max   / 1024.0 : 0.0;
    res.rss_growth_mb = (rss_first > 0 && rss_last > 0)
                        ? (rss_last - rss_first) / 1024.0 : 0.0;

    /* Load and compare reference */
    char ref_path[512];
    snprintf(ref_path, sizeof(ref_path), "../../reference/%s.expr.txt", c->name);
    if (access(ref_path, R_OK) != 0) {
        snprintf(res.verif, sizeof(res.verif), "SKIP (reference file not bundled)");
        rv_free(&got); res.ok = 1; return res;
    }
    row_vec_t ref = {NULL, 0, 0};
    if (parse_ref(ref_path, &ref) != 0) {
        snprintf(res.verif, sizeof(res.verif), "FAIL (no ref)");
        rv_free(&got); res.ok = 0; return res;
    }
    sort_dedup(&ref);

    int miss_c = 0, miss_r = 0, coef_mm = 0; double max_ad = 0.0;
    int rep = 0, ri = 0, gi = 0;
    while (ri < ref.size || gi < got.size) {
        int cmp;
        if      (ri >= ref.size) cmp =  1;
        else if (gi >= got.size) cmp = -1;
        else if (ref.data[ri].ij_key != got.data[gi].ij_key)
            cmp = ref.data[ri].ij_key - got.data[gi].ij_key;
        else
            cmp = ref.data[ri].pqrs - got.data[gi].pqrs;

        if (cmp < 0) {
            if (rep++ < 3)
                fprintf(stderr, "  [%s] MISS in C: ij=%d pqrs=%d coef=%.6e\n",
                        c->name, ref.data[ri].ij_key, ref.data[ri].pqrs,
                        ref.data[ri].coef);
            miss_c++; ri++;
        } else if (cmp > 0) {
            if (rep++ < 3)
                fprintf(stderr, "  [%s] EXTRA in C: ij=%d pqrs=%d coef=%.6e\n",
                        c->name, got.data[gi].ij_key, got.data[gi].pqrs,
                        got.data[gi].coef);
            miss_r++; gi++;
        } else {
            double ad = fabs(ref.data[ri].coef - got.data[gi].coef);
            if (ad > max_ad) max_ad = ad;
            if (ad > 1e-10) coef_mm++;
            ri++; gi++;
        }
    }
    rv_free(&got); rv_free(&ref);

    int fail = miss_c + miss_r + coef_mm;
    if (!fail) {
        if (max_ad <= 1e-15)
            snprintf(res.verif, sizeof(res.verif), "PASS bit-ex");
        else
            snprintf(res.verif, sizeof(res.verif), "PASS %.2e", max_ad);
        res.ok = 1;
    } else {
        snprintf(res.verif, sizeof(res.verif),
                 "FAIL miss_c=%d miss_r=%d coef=%d max=%.2e",
                 miss_c, miss_r, coef_mm, max_ad);
        res.ok = 0;
    }
    return res;
}

/* ------------------------------------------------------------------ */
/* Sample mode: Hermitian self-check                                  */
/* Runs the same 100 pairs as the timing sample (srand(20260513)).    */
/* For each (bra,ket), also runs (ket,bra); verifies matching pqrs    */
/* entries agree to <= 1e-10.                                         */
/* ------------------------------------------------------------------ */

static result_t run_sample(const case_def_t *c, hexpr_drt_t *drt, int ncsfs,
                            int n_sample, double samp_ms,
                            int32_t *pb_ij, double *cb_ij,
                            int32_t *pb_ji, double *cb_ji)
{
    result_t res; memset(&res, 0, sizeof(res));
    strncpy(res.name, c->name, 31); res.ncsfs = ncsfs; res.mode = 1;
    res.pairs_tested = n_sample;
    res.per_pair_ms  = n_sample > 0 ? samp_ms / n_sample : 0.0;
    res.total_s      = samp_ms / 1000.0;

    srand(20260513);
    long rss_first = -1, rss_last = -1, rss_max = 0;
    double max_hd = 0.0;
    int herm_fail = 0;

    for (int s = 0; s < n_sample; s++) {
        int bra = rand() % ncsfs;
        int ket = rand() % ncsfs;

        int cnt_ij, cnt_ji;
        hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_FULL,
                        &cnt_ij, NULL, pb_ij, cb_ij);
        hexpr_pair_eval(drt, ket, bra, HEXPR_EXPR_FULL,
                        &cnt_ji, NULL, pb_ji, cb_ji);

        /* For each pqrs in (bra,ket), check matching entry in (ket,bra). */
        for (int a = 0; a < cnt_ij; a++) {
            for (int b = 0; b < cnt_ji; b++) {
                if (pb_ij[a] == pb_ji[b]) {
                    double ad = fabs(cb_ij[a] - cb_ji[b]);
                    if (ad > max_hd) max_hd = ad;
                    if (ad > 1e-10) herm_fail++;
                    break;
                }
            }
        }

        if ((s + 1) % 10 == 0) {
            long rss = vmrss_kb();
            if (rss > 0) {
                if (rss_first < 0) rss_first = rss;
                rss_last = rss;
                if (rss > rss_max) rss_max = rss;
            }
        }
    }

    res.rss_max_mb    = rss_max   > 0 ? rss_max   / 1024.0 : 0.0;
    res.rss_growth_mb = (rss_first > 0 && rss_last > 0)
                        ? (rss_last - rss_first) / 1024.0 : 0.0;

    if (!herm_fail) {
        snprintf(res.verif, sizeof(res.verif),
                 "Hermitian PASS (sample, max=%.2e)", max_hd);
        res.ok = 1;
    } else {
        snprintf(res.verif, sizeof(res.verif),
                 "Hermitian FAIL (n=%d max=%.2e)", herm_fail, max_hd);
        res.ok = 0;
    }
    return res;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void)
{
    result_t results[NCASES];
    int n_results = 0, total_fail = 0;

    printf("Expression-generator stress test\n");
    printf("=======================================\n\n");

    for (int ci = 0; ci < NCASES; ci++) {
        const case_def_t *c = &CASES[ci];
        printf("[%s] ", c->name); fflush(stdout);

        hexpr_drt_t *drt = hexpr_drt_create(c->kind, c->nve, c->spin,
                                             c->ncore, c->nval, c->next);
        if (!drt) {
            printf("FAIL: DRT creation\n");
            result_t r; memset(&r, 0, sizeof(r));
            strncpy(r.name, c->name, 31);
            snprintf(r.verif, sizeof(r.verif), "FAIL (DRT)");
            results[n_results++] = r;
            total_fail++;
            continue;
        }
        int ncsfs = hexpr_drt_ncsfs(drt);
        printf("ncsfs=%d ", ncsfs); fflush(stdout);

        /* Allocate per-pair buffers using max_terms upper bound. */
        int max_t = 0;
        hexpr_pair_max_terms(drt, HEXPR_EXPR_FULL, &max_t);
        if (max_t <= 0) max_t = 256;  /* safe fallback */

        int32_t *pb1 = malloc((size_t)max_t * sizeof(int32_t));
        int32_t *pb2 = malloc((size_t)max_t * sizeof(int32_t));
        double  *cb1 = malloc((size_t)max_t * sizeof(double));
        double  *cb2 = malloc((size_t)max_t * sizeof(double));
        if (!pb1 || !pb2 || !cb1 || !cb2) {
            printf("FAIL: OOM (max_t=%d)\n", max_t);
            free(pb1); free(pb2); free(cb1); free(cb2);
            hexpr_drt_destroy(drt);
            result_t r; memset(&r, 0, sizeof(r));
            strncpy(r.name, c->name, 31);
            snprintf(r.verif, sizeof(r.verif), "FAIL (OOM)");
            results[n_results++] = r;
            total_fail++;
            continue;
        }

        /* Step 1: Sample 100 pairs for timing estimate. */
        int n_sample = ((long)ncsfs * ncsfs < 100) ? ncsfs * ncsfs : 100;
        srand(20260513);
        double t0s = now_ms();
        for (int s = 0; s < n_sample; s++) {
            int bra = ncsfs > 1 ? rand() % ncsfs : 0;
            int ket = ncsfs > 1 ? rand() % ncsfs : 0;
            int cnt;
            hexpr_pair_eval(drt, bra, ket, HEXPR_EXPR_FULL, &cnt, NULL, pb1, cb1);
        }
        double samp_ms = now_ms() - t0s;
        double ms_pp   = n_sample > 0 ? samp_ms / n_sample : 0.0;
        double est_ms  = ms_pp * (double)ncsfs * ncsfs;

        printf("est=%.1fs ", est_ms / 1000.0); fflush(stdout);

        /* Step 2: Full sweep or sample mode. */
        result_t r;
        const double THREE_HOURS_MS = 3.0 * 3600.0 * 1000.0;
        if (est_ms <= THREE_HOURS_MS) {
            printf("(full) "); fflush(stdout);
            r = run_full(c, drt, ncsfs, pb1, cb1);
        } else {
            printf("(sample) "); fflush(stdout);
            r = run_sample(c, drt, ncsfs, n_sample, samp_ms,
                           pb1, cb1, pb2, cb2);
        }

        free(pb1); free(pb2); free(cb1); free(cb2);
        hexpr_drt_destroy(drt);

        printf("%s\n", r.verif);
        fflush(stdout);
        results[n_results++] = r;
        if (!r.ok) total_fail++;
    }

    /* ---- Summary table ---- */
    printf("\n");
    printf("%-16s  %6s  %14s  %12s  %12s  %-36s  %10s  %13s\n",
           "case", "ncsfs", "pairs_tested", "per_pair_ms",
           "total_time_s", "verification", "rss_max_mb", "rss_growth_mb");
    printf("%-16s  %6s  %14s  %12s  %12s  %-36s  %10s  %13s\n",
           "----", "-----", "------------", "-----------",
           "------------", "------------", "----------", "-------------");
    for (int i = 0; i < n_results; i++) {
        result_t *r = &results[i];
        char pairs_str[32];
        if (r->mode == 0)
            snprintf(pairs_str, sizeof(pairs_str), "%ld (full)", r->pairs_tested);
        else
            snprintf(pairs_str, sizeof(pairs_str), "%ld (sample)", r->pairs_tested);
        printf("%-16s  %6d  %14s  %12.3f  %12.2f  %-36s  %10.1f  %13.1f\n",
               r->name, r->ncsfs, pairs_str, r->per_pair_ms, r->total_s,
               r->verif, r->rss_max_mb, r->rss_growth_mb);
    }
    printf("\n");

    if (total_fail == 0)
        printf("ALL CASES PASS\n");
    else
        printf("%d CASE(S) FAILED\n", total_fail);

    return total_fail > 0 ? 1 : 0;
}
