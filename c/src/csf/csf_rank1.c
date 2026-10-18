/* csf/csf_rank1.c -- reduced matrix elements of one-electron rank-1 spin operators */
/**
 * @file csf_rank1.c
 * @brief @ref hexpr_rank1_pair: <S'||(a+_p a~_q)^1||S> between two CSFs given as
 *        step vectors, by the same per-orbital factors L_j as the spin-free path.
 *
 * Formulation:
 *
 *   R_pq(Psi, Phi) = (-1)^(1 + S' - S) * sqrt(2S + 1) * sigma_pq * prod_j L_j
 *
 * with L_j = hmlkva (expr_tree.c) unchanged, read through
 * hexpr_lk_value_direct. The operator chain, orbitals ascending:
 *
 *   p > q : q carries a~ (OP_A, theta 1/2), p carries a+ (OP_C, theta 1); sigma +1
 *   p < q : p carries a+ (OP_C, theta 1/2), q carries a~ (OP_A, theta 1); sigma -1,
 *           because (a+_p a~_q)^1 = -(a~_q a+_p)^1
 *   p = q : p carries [a+ a~]^1 (OP_CAT, theta 1);                       sigma +1
 *
 * theta is passed to L_j as a multiplicity (2s+1). After the last operator the
 * rank stays 1 up to the top orbital: the walks do not close into a loop.
 * Below the lowest operator the rank is 0 and L_j = 1 when the steps agree,
 * while differing steps there make the element zero. On an orbital without an
 * operator that is empty or doubly occupied in both CSFs L_j = 1 for any rank
 * (the 6j form), so those orbitals are skipped.
 *
 * No DRT is built: a DRT has one head, and the bra and ket may differ in S.
 */
#include "hexpr.h"
#include "expr_internal.h"
#include "hexpr_error.h"

#include <math.h>
#include <stdint.h>

/* Values below this are zero up to rounding of the 9j products. */
#define RANK1_ZERO 1.0e-12

/* Largest 2S the L(k) table key can hold (psi, phi <= 15). */
#define RANK1_MAX_SPIN_X2 14

static const int OCC[4] = { 0, 1, 1, 2 };      /* step e,u,d,f -> electrons */

/* 2S after each orbital (b[j] for orbital j); -1 if the walk is invalid or
 * runs past the key range. */
static int walk_spins(int norb, const unsigned char *s, int *b)
{
    int x = 0;
    for (int j = 0; j < norb; j++) {
        if (s[j] > 3) return -1;
        if (s[j] == 1) x++;
        else if (s[j] == 2) x--;
        if (x < 0 || x > RANK1_MAX_SPIN_X2) return -1;
        b[j] = x;
    }
    return x;
}

/* One reduced element for the chain of (p, q). */
static double eval_pq(int norb, const unsigned char *bra, const unsigned char *ket,
                      const int *bb, const int *bk, int p, int q)
{
    int lo = p < q ? p : q;
    int hi = p < q ? q : p;
    for (int j = 0; j < lo; j++)
        if (bra[j] != ket[j]) return 0.0;      /* rank 0 below: walks must agree */

    double sigma = (p < q) ? -1.0 : 1.0;
    double prod = 1.0;
    int theta = 1;                              /* multiplicity of the chain rank */
    for (int j = lo; j < norb; j++) {
        int hf, th;
        if (p == q && j == p)       { hf = OP_CAT; th = 3; }
        else if (j == lo && p != q) { hf = (p > q) ? OP_A : OP_C; th = 2; }
        else if (j == hi && p != q) { hf = (p > q) ? OP_C : OP_A; th = 3; }
        else {
            hf = OP_UNITY; th = theta;
            if (bra[j] == ket[j] && (bra[j] == 0 || bra[j] == 3))
                continue;                       /* e/e, f/f: L_j = 1 */
        }
        double L = hexpr_lk_value_direct(bra[j] + 1, ket[j] + 1, hf,
                                         bb[j] + 1, bk[j] + 1, th, theta);
        theta = th;
        prod *= L;
        if (prod == 0.0) return 0.0;
    }
    int sb = bb[norb - 1], sk = bk[norb - 1];   /* 2S', 2S */
    double phase = (((2 + sb - sk) / 2) & 1) ? -1.0 : 1.0;   /* (-1)^(1+S'-S) */
    return phase * sqrt((double)(sk + 1)) * sigma * prod;
}

int hexpr_rank1_max_terms(int norb)
{
    return norb > 0 ? norb : 0;
}

hexpr_status_t hexpr_rank1_pair(int norb,
                                const unsigned char *bra_steps,
                                const unsigned char *ket_steps,
                                int *out_count,
                                int32_t *out_p, int32_t *out_q,
                                double *out_R)
{
    static const char *fn = "hexpr_rank1_pair";
    if (out_count) *out_count = 0;
    if (norb <= 0 || bra_steps == NULL || ket_steps == NULL || out_count == NULL) {
        hexpr_set_error("%s: norb must be positive and the step vectors and "
                        "out_count non-NULL", fn);
        return HEXPR_ERR_INVALID_ARG;
    }
    if (hexpr_init() != HEXPR_OK) return HEXPR_ERR_INTERNAL;

    /* 2S after every orbital, for both walks */
    enum { LOCAL = 256 };
    int sbuf[2 * LOCAL];
    int *bb = sbuf, *bk = sbuf + LOCAL;
    if (norb > LOCAL) {
        hexpr_set_error("%s: norb=%d exceeds %d", fn, norb, LOCAL);
        return HEXPR_ERR_RANGE;
    }
    int sb = walk_spins(norb, bra_steps, bb);
    int sk = walk_spins(norb, ket_steps, bk);
    if (sb < 0 || sk < 0) {
        hexpr_set_error("%s: invalid step vector (steps must be 0..3, partial "
                        "2S in 0..%d)", fn, RANK1_MAX_SPIN_X2);
        return HEXPR_ERR_INVALID_ARG;
    }

    /* selection rules: same electron number, |S'-S| <= 1 <= S'+S */
    int nb = 0, nk = 0;
    for (int j = 0; j < norb; j++) { nb += OCC[bra_steps[j]]; nk += OCC[ket_steps[j]]; }
    if (nb != nk || sb - sk > 2 || sk - sb > 2 || sb + sk < 2)
        return HEXPR_OK;

    /* candidates from the occupation difference */
    int p_plus = -1, q_minus = -1, ndiff = 0;
    for (int j = 0; j < norb; j++) {
        int d = OCC[bra_steps[j]] - OCC[ket_steps[j]];
        if (d == 0) continue;
        ndiff++;
        if (d == 1 && p_plus < 0) p_plus = j;
        else if (d == -1 && q_minus < 0) q_minus = j;
        else return HEXPR_OK;                   /* |d| = 2 or a second +1/-1 */
    }

    int n = 0;
    if (ndiff == 2) {
        double R = eval_pq(norb, bra_steps, ket_steps, bb, bk, p_plus, q_minus);
        if (fabs(R) > RANK1_ZERO) {
            if (out_p) out_p[n] = p_plus;
            if (out_q) out_q[n] = q_minus;
            if (out_R) out_R[n] = R;
            n++;
        }
    } else if (ndiff == 0) {
        for (int p = 0; p < norb; p++) {
            if (bra_steps[p] != 1 && bra_steps[p] != 2) continue;   /* open shells */
            double R = eval_pq(norb, bra_steps, ket_steps, bb, bk, p, p);
            if (fabs(R) > RANK1_ZERO) {
                if (out_p) out_p[n] = p;
                if (out_q) out_q[n] = p;
                if (out_R) out_R[n] = R;
                n++;
            }
        }
    }
    *out_count = n;
    return HEXPR_OK;
}
