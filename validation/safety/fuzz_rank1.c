/* libFuzzer harness for hexpr_rank1_pair (libhexpr 1.2): arbitrary bytes as
 * norb and two step vectors (codes outside 0..3, invalid walks and spins
 * included). The call must either return HEXPR_OK with at most
 * hexpr_rank1_max_terms(norb) in-range terms, or an error status; and for
 * valid inputs the adjoint relation R_qp(ket,bra) = (-1)^(S'-S) R_pq(bra,ket)
 * must hold. Run under ASan/UBSan. */
#include <hexpr.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int two_s(const unsigned char *st, int n)
{
    int b = 0;
    for (int i = 0; i < n; i++) b += (st[i] == 1) - (st[i] == 2);
    return b;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    static int inited = 0;
    if (!inited) { hexpr_init(); inited = 1; }
    if (size < 1) return 0;
    int norb = data[0] % 40;                      /* 0..39, beyond the data too */
    unsigned char bra[64] = {0}, ket[64] = {0};
    size_t avail = size - 1;
    for (int i = 0; i < norb; i++) {
        bra[i] = (size_t)i < avail ? data[1 + i] : 0;
        ket[i] = (size_t)(norb + i) < avail ? data[1 + norb + i] : 0;
        /* keep mostly valid codes so the evaluator is reached, sometimes not */
        if (!(data[0] & 0x80)) { bra[i] &= 3; ket[i] &= 3; }
    }
    int nmax = hexpr_rank1_max_terms(norb);
    if (nmax < 1) nmax = 1;
    int32_t *p = malloc(sizeof(int32_t) * nmax), *q = malloc(sizeof(int32_t) * nmax);
    int32_t *p2 = malloc(sizeof(int32_t) * nmax), *q2 = malloc(sizeof(int32_t) * nmax);
    double *r = malloc(sizeof(double) * nmax), *r2 = malloc(sizeof(double) * nmax);
    int n = -1, n2 = -1;
    hexpr_status_t st = hexpr_rank1_pair(norb, bra, ket, &n, p, q, r);
    if (st == HEXPR_OK) {
        if (n < 0 || n > nmax) abort();
        for (int k = 0; k < n; k++)
            if (p[k] < 0 || p[k] >= norb || q[k] < 0 || q[k] >= norb || !isfinite(r[k])) abort();
        if (hexpr_rank1_pair(norb, ket, bra, &n2, p2, q2, r2) != HEXPR_OK || n2 != n) abort();
        int ds = two_s(bra, norb) - two_s(ket, norb);         /* 2(S'-S) */
        double sign = ((ds / 2) % 2 == 0) ? 1.0 : -1.0;
        for (int k = 0; k < n; k++) {
            int found = 0;
            for (int l = 0; l < n2; l++)
                if (p2[l] == q[k] && q2[l] == p[k]) {
                    if (fabs(r2[l] - sign * r[k]) > 1e-12) abort();
                    found = 1;
                }
            if (!found) abort();
        }
    }
    /* NULL output arrays are allowed: count only */
    int n3 = -1;
    if (hexpr_rank1_pair(norb, bra, ket, &n3, NULL, NULL, NULL) != st) abort();
    if (st == HEXPR_OK && n3 != n) abort();
    free(p); free(q); free(p2); free(q2); free(r); free(r2);
    return 0;
}
