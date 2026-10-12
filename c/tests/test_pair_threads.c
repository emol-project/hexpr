/* test_pair_threads.c -- pair evaluations from several threads at once.
 *
 * The evaluation state (the L(k) table and the walk restriction of the v1.1
 * *_full calls) is thread-local, so evaluations may run concurrently: on one
 * shared DRT handle (read-only queries), and on per-thread CSFSET sub-DRTs
 * built inside the workers. For every sampled pair of a few spaces this
 * checks that the terms computed by 8 threads are bit-identical to those of
 * a serial run (same terms, same order, same coefficient bits), for
 * hexpr_pair_eval and hexpr_pair_eval_full. */
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "../include/hexpr.h"

#define NTHREADS 8

typedef struct { int n; int32_t *ij, *pq; double *cf; } res_t;

typedef struct {
    const hexpr_drt_t *drt;
    int npairs;
    const int *bra, *ket;
    int full, sub;              /* *_full call; evaluate on a 2-CSF sub-DRT */
    res_t *out;
    int tid, nth;
    int bad;
} job_t;

static int eval_one(const hexpr_drt_t *drt, int I, int J, int full, int sub, res_t *r)
{
    const hexpr_drt_t *d = drt;
    hexpr_drt_t *own = NULL;
    if (sub) {                  /* as libhexpr2 does: a fresh sub-DRT per pair */
        int norb = hexpr_drt_norb(drt);
        const unsigned char *rows = hexpr_drt_csfs(drt);
        unsigned char buf[2 * 64];
        memcpy(buf, rows + (size_t)I * norb, (size_t)norb);
        memcpy(buf + norb, rows + (size_t)J * norb, (size_t)norb);
        own = hexpr_drt_from_csfset(norb, buf, I == J ? 1 : 2);
        if (own == NULL) return 1;
        d = own;
        J = I == J ? 0 : 1;
        I = 0;
    }
    int mx = 0;
    if (hexpr_pair_max_terms(d, HEXPR_EXPR_FULL, &mx) != HEXPR_OK) goto fail;
    r->ij = malloc((size_t)(mx + 1) * sizeof *r->ij);
    r->pq = malloc((size_t)(mx + 1) * sizeof *r->pq);
    r->cf = malloc((size_t)(mx + 1) * sizeof *r->cf);
    hexpr_status_t st = full
        ? hexpr_pair_eval_full(d, I, J, HEXPR_EXPR_FULL, &r->n, r->ij, r->pq, r->cf)
        : hexpr_pair_eval(d, I, J, HEXPR_EXPR_FULL, &r->n, r->ij, r->pq, r->cf);
    if (own) hexpr_drt_destroy(own);
    return st != HEXPR_OK;
fail:
    if (own) hexpr_drt_destroy(own);
    return 1;
}

static void *worker(void *arg)
{
    job_t *j = arg;
    for (int p = j->tid; p < j->npairs; p += j->nth)
        if (eval_one(j->drt, j->bra[p], j->ket[p], j->full, j->sub, &j->out[p])) j->bad++;
    return NULL;
}

static unsigned rng = 2463534242u;
static int rnd(int n) { rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5; return (int)(rng % (unsigned)n); }

static void free_res(res_t *r, int n)
{
    for (int p = 0; p < n; p++) { free(r[p].ij); free(r[p].pq); free(r[p].cf); }
    free(r);
}

static int run_case(const char *name, const hexpr_drt_t *drt, int npairs)
{
    int n = hexpr_drt_ncsfs(drt), failures = 0;
    int *bra = malloc((size_t)npairs * sizeof *bra), *ket = malloc((size_t)npairs * sizeof *ket);
    for (int p = 0; p < npairs; p++) {   /* half random, half nearby (more terms) */
        bra[p] = rnd(n);
        ket[p] = p % 2 ? rnd(n) : (bra[p] + rnd(8)) % n;
    }
    for (int mode = 0; mode < 4; mode++) {
        int full = mode & 1, sub = mode >> 1;
        res_t *ser = calloc((size_t)npairs, sizeof *ser), *par = calloc((size_t)npairs, sizeof *par);
        job_t s = {drt, npairs, bra, ket, full, sub, ser, 0, 1, 0};
        worker(&s);
        job_t jobs[NTHREADS];
        pthread_t th[NTHREADS];
        for (int t = 0; t < NTHREADS; t++) {
            jobs[t] = (job_t){drt, npairs, bra, ket, full, sub, par, t, NTHREADS, 0};
            pthread_create(&th[t], NULL, worker, &jobs[t]);
        }
        int bad = s.bad;
        for (int t = 0; t < NTHREADS; t++) { pthread_join(th[t], NULL); bad += jobs[t].bad; }
        long nterms = 0;
        for (int p = 0; p < npairs; p++) {
            nterms += ser[p].n;
            if (ser[p].n != par[p].n ||
                memcmp(ser[p].ij, par[p].ij, (size_t)ser[p].n * sizeof(int32_t)) ||
                memcmp(ser[p].pq, par[p].pq, (size_t)ser[p].n * sizeof(int32_t)) ||
                memcmp(ser[p].cf, par[p].cf, (size_t)ser[p].n * sizeof(double))) {
                if (failures < 10)
                    printf("FAIL: %s %s%s pair (%d,%d): threaded result differs\n", name,
                           full ? "eval_full" : "eval", sub ? " on sub-DRT" : "", bra[p], ket[p]);
                failures++;
            }
        }
        if (bad) { printf("FAIL: %s: %d evaluation error(s)\n", name, bad); failures++; }
        printf("  %-14s %-9s %-8s %5d pairs %8ld terms\n", name, full ? "eval_full" : "eval",
               sub ? "sub-DRT" : "shared", npairs, nterms);
        free_res(ser, npairs);
        free_res(par, npairs);
    }
    free(bra); free(ket);
    return failures;
}

int main(void)
{
    if (hexpr_init() != HEXPR_OK) { printf("FAIL: hexpr_init\n"); return 1; }
    printf("=== pair evaluations from %d threads vs serial ===\n", NTHREADS);
    struct { const char *name; hexpr_drt_kind_t k; int nve, s2, nc, nv, ne; } cs[] = {
        {"s_full_6o6e", HEXPR_DRT_SOCI, 6, 0, 0, 6, 0},
        {"t_soci_small", HEXPR_DRT_SOCI, 6, 2, 1, 4, 1},
        {"s_foci_ext2", HEXPR_DRT_FOCI, 6, 0, 1, 3, 3},
    };
    int failures = 0;
    for (size_t c = 0; c < sizeof cs / sizeof cs[0]; c++) {
        hexpr_drt_t *d = hexpr_drt_create(cs[c].k, cs[c].nve, cs[c].s2, cs[c].nc, cs[c].nv, cs[c].ne);
        if (d == NULL) { printf("FAIL: drt %s\n", cs[c].name); failures++; continue; }
        failures += run_case(cs[c].name, d, 2000);
        hexpr_drt_destroy(d);
    }
    if (failures) { printf("FAILED: %d check(s)\n", failures); return 1; }
    printf("ALL PASS\n");
    return 0;
}
