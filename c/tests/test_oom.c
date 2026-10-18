/* test_oom.c -- allocation-failure injection over the public API (since 1.2).
 *
 * Linked against libhexpr.a with -Wl,--wrap=malloc,--wrap=calloc,--wrap=realloc
 * (GNU ld), so every allocation the library makes goes through the wrappers
 * below. For each scenario the test first runs it normally, records a hash of
 * its results and the number of allocations it makes (A), and then reruns it
 * A+1 times with the k-th allocation failing (k = 1..A+1). Each faulted run
 * must either
 *   - report out of memory: HEXPR_ERR_OOM from a status-returning call, or
 *     NULL with "out of memory" in hexpr_last_error() from a builder; or
 *   - succeed with exactly the results of the normal run;
 * and the normal run repeated right after must still give those results (no
 * state left behind). Crashes, other statuses and wrong results fail the
 * test; under ASan/LSan (validation/safety/run_matrix.sh) so do invalid
 * accesses and leaks on the failure paths.
 *
 * The test itself allocates nothing while a fault is armed (static buffers). */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/hexpr.h"

/* ---- allocation wrappers ------------------------------------------------ */
void *__real_malloc(size_t n);
void *__real_calloc(size_t n, size_t s);
void *__real_realloc(void *p, size_t n);

static long g_count;          /* allocations seen since the last arm/reset */
static long g_fail_at = -1;   /* fail this allocation (1-based); -1 = never */

static int fail_now(void) { return ++g_count == g_fail_at; }
void *__wrap_malloc(size_t n)             { return fail_now() ? NULL : __real_malloc(n); }
void *__wrap_calloc(size_t n, size_t s)   { return fail_now() ? NULL : __real_calloc(n, s); }
void *__wrap_realloc(void *p, size_t n)   { return fail_now() ? NULL : __real_realloc(p, n); }

/* ---- result hashing ----------------------------------------------------- */
static uint64_t h_bytes(uint64_t h, const void *p, size_t n)
{
    const unsigned char *b = p;
    for (size_t i = 0; i < n; i++) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}
#define H0 14695981039346656037ULL
#define HV(h, v) h_bytes((h), &(v), sizeof(v))

/* A scenario returns 0 with *out_hash set, 1 if it reported out of memory
 * correctly, or -1 for anything else (wrong status, NULL without the OOM
 * message). */
typedef int (*scenario_fn)(uint64_t *out_hash);

static int is_oom_msg(void) { return strstr(hexpr_last_error(), "out of memory") != NULL; }

static int status_result(hexpr_status_t st)
{
    if (st == HEXPR_OK) return 0;
    if (st == HEXPR_ERR_OOM) return 1;
    fprintf(stderr, "    unexpected status %d: %s\n", (int)st, hexpr_last_error());
    return -1;
}

/* Fixtures, built once before any fault is armed. */
static hexpr_drt_t *g_soci, *g_foci;

static uint64_t hash_drt(const hexpr_drt_t *d)
{
    uint64_t h = H0;
    int n = hexpr_drt_ncsfs(d), norb = hexpr_drt_norb(d), nn = hexpr_drt_nnodes(d);
    int nw = hexpr_drt_walk_ncsfs(d);
    h = HV(h, n); h = HV(h, norb); h = HV(h, nn); h = HV(h, nw);
    h = h_bytes(h, hexpr_drt_csfs(d), (size_t)n * (size_t)norb);
    h = h_bytes(h, hexpr_drt_walk_csfs(d), (size_t)nw * (size_t)norb);
    return h;
}

static int sc_create_soci(uint64_t *oh)
{
    hexpr_drt_t *d = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 3, 2);
    if (!d) return is_oom_msg() ? 1 : -1;
    *oh = hash_drt(d);
    hexpr_drt_destroy(d);
    return 0;
}

static int sc_create_foci(uint64_t *oh)
{
    hexpr_drt_t *d = hexpr_drt_create(HEXPR_DRT_FOCI, 5, 1, 1, 2, 2);
    if (!d) return is_oom_msg() ? 1 : -1;
    *oh = hash_drt(d);
    hexpr_drt_destroy(d);
    return 0;
}

static unsigned char g_rows[3 * 64];
static int sc_from_csfset(uint64_t *oh)
{
    int norb = hexpr_drt_norb(g_soci);
    hexpr_drt_t *d = hexpr_drt_from_csfset(norb, g_rows, 3);
    if (!d) return is_oom_msg() ? 1 : -1;
    *oh = hash_drt(d);
    hexpr_drt_destroy(d);
    return 0;
}

static int32_t g_ij[4096], g_pq[4096], g_bp[4096], g_kp[4096];
static double  g_cf[4096];

static int sc_expr_build(uint64_t *oh)
{
    hexpr_expr_t *e = hexpr_expr_build(g_foci);
    if (!e) return is_oom_msg() ? 1 : -1;
    int64_t n = hexpr_expr_nterms_full(e);
    uint64_t h = H0;
    h = HV(h, n);
    if (n > 4096 || hexpr_expr_read(e, HEXPR_EXPR_FULL, g_ij, g_pq, g_cf) != HEXPR_OK) {
        hexpr_expr_destroy(e);
        return -1;
    }
    h = h_bytes(h, g_ij, (size_t)n * sizeof *g_ij);
    h = h_bytes(h, g_pq, (size_t)n * sizeof *g_pq);
    h = h_bytes(h, g_cf, (size_t)n * sizeof *g_cf);
    hexpr_expr_destroy(e);
    *oh = h;
    return 0;
}

static uint64_t hash_terms(int n, int with_pos)
{
    uint64_t h = H0;
    h = HV(h, n);
    h = h_bytes(h, g_ij, (size_t)n * sizeof *g_ij);
    h = h_bytes(h, g_pq, (size_t)n * sizeof *g_pq);
    h = h_bytes(h, g_cf, (size_t)n * sizeof *g_cf);
    if (with_pos) {
        h = h_bytes(h, g_bp, (size_t)n * sizeof *g_bp);
        h = h_bytes(h, g_kp, (size_t)n * sizeof *g_kp);
    }
    return h;
}

static int sc_pair_eval(uint64_t *oh)
{
    int n = -1;
    int r = status_result(hexpr_pair_eval(g_soci, 3, 7, HEXPR_EXPR_FULL, &n, g_ij, g_pq, g_cf));
    if (r == 0) *oh = hash_terms(n, 0);
    return r;
}

static int sc_pair_eval_full(uint64_t *oh)
{
    int n = -1;
    int r = status_result(hexpr_pair_eval_full(g_soci, 7, 3, HEXPR_EXPR_FULL, &n, g_ij, g_pq, g_cf));
    if (r == 0) *oh = hash_terms(n, 0);
    return r;
}

static int sc_pair_count(uint64_t *oh)
{
    int n = -1;
    int r = status_result(hexpr_pair_count_full(g_soci, 2, 5, HEXPR_EXPR_ONE, &n));
    if (r == 0) { uint64_t h = H0; *oh = HV(h, n); }
    return r;
}

static int sc_block_full(uint64_t *oh)
{
    static const int bra[2] = {0, 4}, ket[3] = {1, 4, 6};
    int n = -1;
    int r = status_result(hexpr_block_eval_full(g_soci, bra, 2, ket, 3, HEXPR_EXPR_FULL, &n,
                                                g_bp, g_kp, g_ij, g_pq, g_cf));
    if (r == 0) *oh = hash_terms(n, 1);
    return r;
}

static int sc_subspace(uint64_t *oh)
{
    static const int idx[3] = {1, 2, 8};
    int n = -1;
    int r = status_result(hexpr_subspace_eval(g_foci, idx, 3, HEXPR_EXPR_FULL, &n,
                                              g_bp, g_kp, g_ij, g_pq, g_cf));
    if (r == 0) *oh = hash_terms(n, 1);
    return r;
}

static uint64_t g_al[256], g_be[256];
static int sc_csf_to_dets(uint64_t *oh)
{
    int n = -1;
    int r = status_result(hexpr_csf_to_dets(g_foci, 5, 1, &n, g_al, g_be, g_cf));
    if (r == 0) {
        uint64_t h = H0;
        h = HV(h, n);
        h = h_bytes(h, g_al, (size_t)n * sizeof *g_al);
        h = h_bytes(h, g_be, (size_t)n * sizeof *g_be);
        *oh = h_bytes(h, g_cf, (size_t)n * sizeof *g_cf);
    }
    return r;
}

/* ---- driver ------------------------------------------------------------- */
static int failures = 0;

static void run(const char *name, scenario_fn fn)
{
    uint64_t h0 = 0, h = 0;
    g_fail_at = -1; g_count = 0;
    if (fn(&h0) != 0) {
        printf("  %-22s FAIL (the normal run did not succeed)\n", name);
        failures++;
        return;
    }
    long nalloc = g_count;
    long oom = 0, ok = 0, bad = 0;
    for (long k = 1; k <= nalloc + 1; k++) {
        g_count = 0; g_fail_at = k;
        int r = fn(&h);
        g_fail_at = -1;
        if (r == 1) oom++;
        else if (r == 0 && h == h0) ok++;
        else {
            bad++;
            if (bad <= 3)
                printf("    %s: allocation %ld -> %s\n", name, k,
                       r == 0 ? "success with different results" : "wrong error report");
        }
        uint64_t h1 = 0;
        g_count = 0;
        if (fn(&h1) != 0 || h1 != h0) {
            bad++;
            if (bad <= 3)
                printf("    %s: the normal run after failing allocation %ld differs\n",
                       name, k);
        }
    }
    printf("  %-22s %6ld allocations: %6ld out of memory, %6ld unaffected  %s\n",
           name, nalloc, oom, ok, bad ? "FAIL" : "PASS");
    if (bad) failures++;
}

int main(void)
{
    hexpr_init();
    g_soci = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 3, 2);
    g_foci = hexpr_drt_create(HEXPR_DRT_FOCI, 5, 1, 1, 2, 2);
    if (!g_soci || !g_foci || hexpr_drt_ncsfs(g_soci) < 9 || hexpr_drt_ncsfs(g_foci) < 9) {
        printf("test_oom: fixtures could not be built\n");
        return 1;
    }
    int norb = hexpr_drt_norb(g_soci);
    for (int r = 0; r < 3; r++)
        memcpy(g_rows + r * norb, hexpr_drt_csfs(g_soci) + (size_t)(2 * r + 1) * norb,
               (size_t)norb);

    printf("=== test_oom: failing each allocation in turn ===\n");
    run("drt_create SOCI",   sc_create_soci);
    run("drt_create FOCI",   sc_create_foci);
    run("drt_from_csfset",   sc_from_csfset);
    run("expr_build",        sc_expr_build);
    run("pair_eval",         sc_pair_eval);
    run("pair_eval_full",    sc_pair_eval_full);
    run("pair_count_full",   sc_pair_count);
    run("block_eval_full",   sc_block_full);
    run("subspace_eval",     sc_subspace);
    run("csf_to_dets",       sc_csf_to_dets);

    hexpr_drt_destroy(g_soci);
    hexpr_drt_destroy(g_foci);
    printf("test_oom: %s\n", failures ? "FAIL" : "PASS");
    return failures ? 1 : 0;
}
