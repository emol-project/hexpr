/* test_rank1.c -- hexpr_rank1_pair, checks that need no external program.
 *
 * The values themselves were settled against determinants first (during
 * development: every pair and (p,q) up to 4 orbitals, random pairs at 6 and
 * 8 orbitals, all within 5e-15). This file keeps what can be checked from
 * inside:
 *
 *   1. total spin: sum_p R(I,I,p,p) = -sqrt2 sqrt(S(S+1)(2S+1)) for every CSF,
 *      i.e. sum_p (a+_p a~_p)^1 = -sqrt2 S;
 *   2. adjoint: R_qp(ket,bra) = (-1)^(S'-S) R_pq(bra,ket) for every pair,
 *      including pairs of different S;
 *   3. worked examples of the per-orbital factors and a singlet-triplet
 *      spin density pair;
 *   4. selection rules give no terms; bad arguments are rejected;
 *   5. at most hexpr_rank1_max_terms(norb) terms.
 */
#include "hexpr.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NORB_MAX 6
#define TOL 1e-12

static int failures = 0;

#define CHECK(cond, ...) do { if (!(cond)) { failures++; \
    fprintf(stderr, "FAIL %s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); } } while (0)

/* all step vectors with nel electrons, any final spin */
typedef struct { unsigned char s[NORB_MAX]; int s2; } csf_t;

static int enumerate(int norb, int nel, csf_t *out, int cap)
{
    int n = 0, total = 1;
    for (int j = 0; j < norb; j++) total *= 4;
    for (int code = 0; code < total; code++) {
        csf_t c; int x = code, e = 0, b = 0, ok = 1;
        for (int j = 0; j < norb; j++) {
            c.s[j] = (unsigned char)(x & 3); x >>= 2;
            e += (c.s[j] == 0) ? 0 : (c.s[j] == 3 ? 2 : 1);
            b += (c.s[j] == 1) ? 1 : (c.s[j] == 2 ? -1 : 0);
            if (b < 0) { ok = 0; break; }
        }
        if (!ok || e != nel) continue;
        c.s2 = b;
        if (n < cap) out[n] = c;
        n++;
    }
    return n;
}

/* R_pq for one pair, as a dense norb x norb table */
static void dense(int norb, const unsigned char *bra, const unsigned char *ket,
                  double *R)
{
    int32_t p[NORB_MAX], q[NORB_MAX];
    double r[NORB_MAX];
    int n = -1;
    memset(R, 0, sizeof(double) * norb * norb);
    hexpr_status_t st = hexpr_rank1_pair(norb, bra, ket, &n, p, q, r);
    CHECK(st == HEXPR_OK, "status %d: %s", (int)st, hexpr_last_error());
    CHECK(n >= 0 && n <= hexpr_rank1_max_terms(norb), "count %d", n);
    for (int k = 0; k < n; k++) R[p[k] * norb + q[k]] = r[k];
}

static double one(const char *bra, const char *ket, int p, int q)
{
    int norb = (int)strlen(bra);
    unsigned char b[NORB_MAX], k[NORB_MAX];
    for (int j = 0; j < norb; j++) {
        b[j] = (unsigned char)(strchr("eudf", bra[j]) - "eudf");
        k[j] = (unsigned char)(strchr("eudf", ket[j]) - "eudf");
    }
    double R[NORB_MAX * NORB_MAX];
    dense(norb, b, k, R);
    return R[p * norb + q];
}

static void identities(int norb)
{
    enum { CAP = 4096 };
    csf_t *c = malloc(sizeof(csf_t) * CAP);
    double *R1 = malloc(sizeof(double) * NORB_MAX * NORB_MAX);
    double *R2 = malloc(sizeof(double) * NORB_MAX * NORB_MAX);
    double wspin = 0.0, wadj = 0.0;
    long npairs = 0;
    for (int nel = 1; nel < 2 * norb; nel++) {
        int n = enumerate(norb, nel, c, CAP);
        CHECK(n <= CAP, "enumeration overflow");
        for (int i = 0; i < n; i++) {
            if (c[i].s2 > 0) {
                dense(norb, c[i].s, c[i].s, R1);
                double tot = 0.0, s = c[i].s2 / 2.0;
                for (int p = 0; p < norb; p++) tot += R1[p * norb + p];
                double ref = -sqrt(2.0) * sqrt(s * (s + 1) * (2 * s + 1));
                if (fabs(tot - ref) > wspin) wspin = fabs(tot - ref);
            }
            for (int k = 0; k < n; k++) {
                int d = c[i].s2 - c[k].s2;
                if (d > 2 || d < -2) continue;
                dense(norb, c[i].s, c[k].s, R1);
                dense(norb, c[k].s, c[i].s, R2);
                double sg = ((d / 2) & 1) ? -1.0 : 1.0;    /* (-1)^(S'-S) */
                for (int p = 0; p < norb; p++)
                    for (int q = 0; q < norb; q++) {
                        double e = fabs(R2[q * norb + p] - sg * R1[p * norb + q]);
                        if (e > wadj) wadj = e;
                    }
                npairs++;
            }
        }
    }
    printf("  norb=%d: %ld pairs, total spin %.1e, adjoint %.1e\n",
           norb, npairs, wspin, wadj);
    CHECK(wspin < TOL, "total-spin identity, norb=%d: %.3e", norb, wspin);
    CHECK(wadj < TOL, "adjoint relation, norb=%d: %.3e", norb, wadj);
    free(c); free(R1); free(R2);
}

int main(void)
{
    CHECK(hexpr_init() == HEXPR_OK, "hexpr_init");

    /* 1, 2 */
    for (int norb = 2; norb <= 5; norb++) identities(norb);

    /* 3: worked examples (bra, ket, p, q -> R), 0-based orbitals */
    struct { const char *bra, *ket; int p, q; double R; } ex[] = {
        { "udf", "fuu", 2, 0,  sqrt(1.5) },   /* Table spin-ex-A: (a+_3 a~_1)^1 */
        { "uuu", "udu", 0, 0,  sqrt(2.0) },   /* Table spin-ex-C: [a+_1 a~_1]^1 */
        { "ud",  "uu",  0, 0, -sqrt(1.5) },   /* singlet <- triplet spin density */
        { "ud",  "uu",  1, 1,  sqrt(1.5) },   /*   ... sums to zero across S     */
        { "uu",  "ud",  0, 0,  sqrt(1.5) },   /* reverse of row 3: (-1)^(0-1) */
    };
    for (size_t i = 0; i < sizeof ex / sizeof ex[0]; i++) {
        double v = one(ex[i].bra, ex[i].ket, ex[i].p, ex[i].q);
        CHECK(fabs(v - ex[i].R) < TOL, "example %zu: %s %s (%d,%d) = %.15f, want %.15f",
              i, ex[i].bra, ex[i].ket, ex[i].p, ex[i].q, v, ex[i].R);
    }

    /* 4: selection rules and argument checks */
    {
        unsigned char a[3] = {1, 2, 0}, b[3] = {0, 1, 2};   /* ud e / e ud: two singlets */
        int n = -1;
        CHECK(hexpr_rank1_pair(3, a, b, &n, NULL, NULL, NULL) == HEXPR_OK && n == 0,
              "two singlets must give no terms (n=%d)", n);
        unsigned char c[2] = {3, 0}, d[2] = {0, 3};         /* fe / ef: two electrons moved */
        CHECK(hexpr_rank1_pair(2, c, d, &n, NULL, NULL, NULL) == HEXPR_OK && n == 0,
              "two-electron move must give no terms (n=%d)", n);
        unsigned char e[2] = {1, 1}, f[2] = {1, 0};         /* different electron number */
        CHECK(hexpr_rank1_pair(2, e, f, &n, NULL, NULL, NULL) == HEXPR_OK && n == 0,
              "different N must give no terms (n=%d)", n);
        unsigned char g[2] = {2, 1}, h[2] = {4, 0};         /* d first: invalid; step 4 */
        CHECK(hexpr_rank1_pair(2, g, e, &n, NULL, NULL, NULL) == HEXPR_ERR_INVALID_ARG,
              "negative partial spin must be rejected");
        CHECK(hexpr_rank1_pair(2, h, e, &n, NULL, NULL, NULL) == HEXPR_ERR_INVALID_ARG,
              "step 4 must be rejected");
        CHECK(hexpr_rank1_pair(0, e, e, &n, NULL, NULL, NULL) == HEXPR_ERR_INVALID_ARG,
              "norb 0 must be rejected");
        CHECK(hexpr_rank1_pair(2, e, e, NULL, NULL, NULL, NULL) == HEXPR_ERR_INVALID_ARG,
              "NULL out_count must be rejected");
    }

    if (failures == 0) printf("test_rank1: PASS\n");
    else printf("test_rank1: %d FAILURE(S)\n", failures);
    return failures ? 1 : 0;
}
