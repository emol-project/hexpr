/* test_csf2det.c -- CSF -> determinant expansion.
 *
 * Two stages, written in that order on purpose.
 *
 * STAGE ONE, the checks that need no oracle: determinant count, coefficient
 * magnitudes, unit norm, Sz, <I|J> = delta_IJ across every Ms, and the
 * orbital direction the bit mapping assumes. Every one of these is a
 * quadratic form or a bit count, so every one of them is blind to the sign
 * of a coefficient. They all passed while csf2det.c was coupling in the
 * wrong order -- see the Eq. (5) note in that file -- which is exactly the
 * limit of what a self-check can do.
 *
 * STAGE TWO, added only after validation/csf2det_vs_pyscf.py had settled the
 * convention against an independent determinant-basis code:
 *
 *   - S^2 in the determinant basis must return S(S+1) on every expanded CSF.
 *     This is sensitive to the coupling order in a way stage one is not,
 *     because S_+ and S_- move one electron between spin blocks and pick up
 *     an operator-ordering sign while doing it.
 *
 *   - the relative sign of the two determinants of an open-shell singlet and
 *     of the corresponding Ms=0 triplet, which stage one could only observe
 *     and print.
 *
 * Writing S^2 here was deliberately deferred. A determinant-basis S^2 repeats
 * the same operator-ordering reasoning csf2det.c does, so writing it first
 * would have risked two implementations sharing one mistake and agreeing.
 * The deferral paid off: two successive hand derivations of the S_+ sign were
 * written here and both were wrong -- the first made S_+ |psi> cancel to
 * zero, the second got the magnitude right and the sign backwards -- and
 * because the expected numbers were already fixed by PySCF, each failure was
 * attributable to THIS file rather than to csf2det.c.
 *
 * The sign was then settled by enumeration rather than by a third derivation.
 * See the note above spin_move_sign().
 *
 * MULTI-WORD CASES. csf2det.c packs determinants into (norb + 63) / 64 words.
 * Everything else in this file uses DRTs with norb <= 5, i.e. one word, so
 * the packing was at one point written and never executed. The DRTs named
 * "wide" below have norb = 69 or 70 and exercise two words, including the
 * partly-filled second word where bit p of word w is orbital 64*w + p.
 *
 * Nothing here hardcodes an expansion that only this library could produce;
 * each check is a hand value, an internal consistency relation, or a
 * quantity whose value was fixed by the PySCF comparison and is cited as
 * such.
 */
#include "hexpr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TOL 1e-12

/* Upper bound on the words this file's S^2 machinery handles. The
 * IMPLEMENTATION is general; only the test is capped, because carrying a
 * variable-length determinant through the S_- S_+ term list would make this
 * file materially harder to read, and this file has already been the source
 * of two wrong signs. Two words are enough: the bit arithmetic p >> 6 and
 * p & 63 does not depend on which word p lands in, so a second word that is
 * only partly filled exercises exactly the same code path a third would. */
#define MAXW 2

static int failures = 0;

static void ok(int cond, const char *what)
{
    printf("  %-58s %s\n", what, cond ? "PASS" : "FAIL");
    if (!cond) failures++;
}

static void ok_close(double got, double want, const char *what)
{
    int cond = fabs(got - want) < 1e-10;
    printf("  %-58s %s", what, cond ? "PASS" : "FAIL");
    if (!cond) { printf("   (got %.15g want %.15g)", got, want); failures++; }
    putchar('\n');
}

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

static int popcount_words(const uint64_t *w, int nwords)
{
    int n = 0;
    for (int i = 0; i < nwords; i++) {
        uint64_t v = w[i];
        while (v) { v &= v - 1; n++; }
    }
    return n;
}

static int words_equal(const uint64_t *a, const uint64_t *b, int nwords)
{
    for (int i = 0; i < nwords; i++) if (a[i] != b[i]) return 0;
    return 1;
}

static int bit_of(const uint64_t *w, int p)
{
    return (int)((w[p >> 6] >> (p & 63)) & 1u);
}

static void bit_set(uint64_t *w, int p)   { w[p >> 6] |=  ((uint64_t)1 << (p & 63)); }
static void bit_clear(uint64_t *w, int p) { w[p >> 6] &= ~((uint64_t)1 << (p & 63)); }

/* Expand one CSF; caller frees pa, pb and pc. Returns ndets, or -1. */
static int expand(const hexpr_drt_t *drt, int index, int ms_x2,
                  uint64_t **pa, uint64_t **pb, double **pc)
{
    int ndets = 0;
    if (hexpr_csf_ndets(drt, index, ms_x2, &ndets) != HEXPR_OK) return -1;

    const int nwords = (hexpr_drt_norb(drt) + 63) / 64;
    uint64_t *a = malloc(sizeof(uint64_t) * (size_t)ndets * (size_t)nwords);
    uint64_t *b = malloc(sizeof(uint64_t) * (size_t)ndets * (size_t)nwords);
    double   *c = malloc(sizeof(double)   * (size_t)ndets);
    if (!a || !b || !c) { free(a); free(b); free(c); return -1; }

    int got = 0;
    if (hexpr_csf_to_dets(drt, index, ms_x2, &got, a, b, c) != HEXPR_OK ||
        got != ndets) {
        free(a); free(b); free(c);
        return -1;
    }
    *pa = a; *pb = b; *pc = c;
    return ndets;
}

/* Print an expansion as orbital-occupation strings, so a small case can be
 * checked by eye. Code per orbital: '.' empty, 'a' alpha, 'b' beta, '2' both. */
static void print_expansion(const char *label, const hexpr_drt_t *drt,
                            const uint64_t *a, const uint64_t *b,
                            const double *c, int ndets)
{
    const int norb = hexpr_drt_norb(drt);
    const int nw   = (norb + 63) / 64;
    printf("     %s:\n", label);
    for (int k = 0; k < ndets; k++) {
        printf("       det%-2d ", k);
        for (int p = 0; p < norb; p++) {
            int ia = bit_of(a + (size_t)k * nw, p);
            int ib = bit_of(b + (size_t)k * nw, p);
            putchar(ia && ib ? '2' : ia ? 'a' : ib ? 'b' : '.');
        }
        printf("   c = %+.15f\n", c[k]);
    }
}

/* <I|J> through the determinant basis. Determinants are orthonormal, so the
 * overlap is the sum of c_I(k) c_J(l) over coinciding (alpha,beta) pairs. */
static double overlap(const hexpr_drt_t *drt, int i, int j, int ms_x2, int *bad)
{
    uint64_t *ai, *bi, *aj, *bj; double *ci, *cj;
    int ni = expand(drt, i, ms_x2, &ai, &bi, &ci);
    if (ni < 0) { *bad = 1; return 0.0; }
    int nj = expand(drt, j, ms_x2, &aj, &bj, &cj);
    if (nj < 0) { free(ai); free(bi); free(ci); *bad = 1; return 0.0; }

    const int nw = (hexpr_drt_norb(drt) + 63) / 64;
    double s = 0.0;
    for (int k = 0; k < ni; k++) {
        if (ci[k] == 0.0) continue;
        for (int l = 0; l < nj; l++) {
            if (cj[l] == 0.0) continue;
            if (words_equal(ai + (size_t)k * nw, aj + (size_t)l * nw, nw) &&
                words_equal(bi + (size_t)k * nw, bj + (size_t)l * nw, nw))
                s += ci[k] * cj[l];
        }
    }
    free(ai); free(bi); free(ci);
    free(aj); free(bj); free(cj);
    return s;
}

/* Find the CSF of `drt` whose step vector matches `want` ("udee" etc.).
 * Returns its index, or -1. */
static int find_csf(const hexpr_drt_t *drt, const char *want)
{
    const int norb  = hexpr_drt_norb(drt);
    const int ncsfs = hexpr_drt_ncsfs(drt);
    static const char CH[4] = { 'e', 'u', 'd', 'f' };
    unsigned char *step = malloc((size_t)norb);
    if (!step) return -1;
    if ((int)strlen(want) != norb) { free(step); return -1; }

    for (int i = 0; i < ncsfs; i++) {
        if (hexpr_csf_steps(drt, i, step) != HEXPR_OK) continue;
        int same = 1;
        for (int k = 0; k < norb; k++)
            if (CH[step[k]] != want[k]) { same = 0; break; }
        if (same) { free(step); return i; }
    }
    free(step);
    return -1;
}

/* ================================================================== */
/* STAGE ONE                                                           */
/* ================================================================== */

/* ------------------------------------------------------------------ */
/* (1) orbital direction                                               */
/* ------------------------------------------------------------------ */
/* With ncore=1, nval=1, next=2 and norb=4, SOCI admits a step vector with
 * both electrons in the LAST position (two external electrons) while FOCI
 * does not; FOCI instead admits both in the FIRST position. If the step
 * vector ran the other way round, external would be at the front and the two
 * kinds would swap. csf2det.c's orbital mapping depends on this. */
static int has_steps(const hexpr_drt_t *drt, const char *want)
{
    return find_csf(drt, want) >= 0;
}

static void test_orbital_direction(void)
{
    printf("-- (1) orbital direction: step[0] is the FIRST orbital --\n");

    hexpr_drt_t *soci = hexpr_drt_create(HEXPR_DRT_SOCI, 2, 0, 1, 1, 2);
    hexpr_drt_t *foci = hexpr_drt_create(HEXPR_DRT_FOCI, 2, 0, 1, 1, 2);
    if (!soci || !foci) {
        printf("  FAIL: hexpr_drt_create returned NULL (%s)\n", hexpr_last_error());
        failures++;
        if (soci) hexpr_drt_destroy(soci);
        if (foci) hexpr_drt_destroy(foci);
        return;
    }

    ok(has_steps(soci, "eeef"),  "SOCI admits eeef (two electrons, last orbital)");
    ok(!has_steps(foci, "eeef"), "FOCI rejects eeef (would be two external)");
    ok(has_steps(foci, "feee"),  "FOCI admits feee (two electrons, first orbital)");

    hexpr_drt_destroy(soci);
    hexpr_drt_destroy(foci);
}

/* ------------------------------------------------------------------ */
/* (2)(5) the two-open-shell cases, magnitudes and relative signs       */
/* ------------------------------------------------------------------ */

/* Two open shells coupled to S=0, expanded at Ms=0.
 *
 * Magnitudes and norm follow from the step vector alone. The RELATIVE sign
 * of the two coefficients does not: it is fixed by the operator ordering,
 * and it was settled against PySCF (validation/csf2det_vs_pyscf.py, H2
 * CAS(2,2) and H4 CAS(4,4), matching to 2e-15). The value recorded there is
 * SAME sign for the singlet.
 *
 * Note what is NOT asserted: the coefficients come out both negative here,
 * but an overall sign common to every determinant of a CSF is a phase that
 * no comparison in the CSF or determinant basis can fix -- <I|J> and
 * U H_det U^T are both invariant under it. Asserting it would pin an
 * arbitrary choice as though it had been measured. */
static void test_singlet_pair(void)
{
    printf("-- (2) CG closed forms, singlet of two open shells --\n");

    hexpr_drt_t *drt = hexpr_drt_create(HEXPR_DRT_FOCI, 2, 0, 1, 1, 2);
    if (!drt) { printf("  FAIL: create (%s)\n", hexpr_last_error()); failures++; return; }

    int idx = find_csf(drt, "udee");
    if (idx < 0) {
        printf("  FAIL: no udee CSF found\n"); failures++;
        hexpr_drt_destroy(drt); return;
    }

    uint64_t *a, *b; double *c;
    int nd = expand(drt, idx, 0, &a, &b, &c);
    ok(nd == 2, "udee at Ms=0 expands into exactly 2 determinants");
    if (nd == 2) {
        ok_close(fabs(c[0]), 1.0 / sqrt(2.0), "|c[0]| == 1/sqrt(2)");
        ok_close(fabs(c[1]), 1.0 / sqrt(2.0), "|c[1]| == 1/sqrt(2)");
        ok_close(c[0] * c[0] + c[1] * c[1], 1.0, "normalized");
        print_expansion("udee S=0 Ms=0", drt, a, b, c, nd);
        ok(c[0] * c[1] > 0.0,
           "singlet: SAME relative sign (fixed by the PySCF comparison)");
    }
    if (nd > 0) { free(a); free(b); free(c); }
    hexpr_drt_destroy(drt);
}

/* The same two orbitals coupled to S=1, expanded at Ms=0.
 *
 * Two assertions of different standing. That the singlet's and the triplet's
 * relative signs are opposite to each other follows from orthogonality --
 * they span the same two determinants -- and holds under ANY operator
 * ordering. That the triplet's own relative sign is OPPOSITE is the value
 * fixed by the PySCF comparison. */
static void test_triplet_pair(void)
{
    printf("-- (5) singlet/triplet relative signs at Ms=0 --\n");

    hexpr_drt_t *sing = hexpr_drt_create(HEXPR_DRT_FOCI, 2, 0, 1, 1, 2);
    hexpr_drt_t *trip = hexpr_drt_create(HEXPR_DRT_FOCI, 2, 2, 1, 1, 2);
    if (!sing || !trip) {
        printf("  FAIL: create (%s)\n", hexpr_last_error()); failures++;
        if (sing) hexpr_drt_destroy(sing);
        if (trip) hexpr_drt_destroy(trip);
        return;
    }

    int is_ = find_csf(sing, "udee");
    int it  = find_csf(trip, "uuee");
    if (is_ < 0 || it < 0) {
        printf("  SKIP: udee or uuee not present in these DRTs\n");
        hexpr_drt_destroy(sing); hexpr_drt_destroy(trip);
        return;
    }

    uint64_t *as, *bs, *at, *bt; double *cs, *ct;
    int ns = expand(sing, is_, 0, &as, &bs, &cs);
    int nt = expand(trip, it,  0, &at, &bt, &ct);

    ok(ns == 2 && nt == 2, "udee and uuee each expand into 2 determinants at Ms=0");
    if (ns == 2 && nt == 2) {
        ok_close(fabs(ct[0]), 1.0 / sqrt(2.0), "|triplet c[0]| == 1/sqrt(2)");
        ok_close(fabs(ct[1]), 1.0 / sqrt(2.0), "|triplet c[1]| == 1/sqrt(2)");

        print_expansion("udee S=0 Ms=0", sing, as, bs, cs, ns);
        print_expansion("uuee S=1 Ms=0", trip, at, bt, ct, nt);

        ok((cs[0] * cs[1]) * (ct[0] * ct[1]) < 0.0,
           "singlet and triplet relative signs are opposite (orthogonality)");
        ok(ct[0] * ct[1] < 0.0,
           "triplet: OPPOSITE relative sign (fixed by the PySCF comparison)");
    }
    if (ns > 0) { free(as); free(bs); free(cs); }
    if (nt > 0) { free(at); free(bt); free(ct); }
    hexpr_drt_destroy(sing);
    hexpr_drt_destroy(trip);
}

/* ================================================================== */
/* STAGE TWO -- S^2 in the determinant basis                           */
/* ================================================================== */

/* Sign picked up when one electron changes spin at orbital p, in the
 * all-alpha-then-all-beta ordering:
 *
 *     sign = (-1)^(alphas below p + betas below p)
 *
 * i.e. the parity of the number of electrons in orbitals strictly below p,
 * counting both spins. The same expression serves S_+ and S_-.
 *
 * HOW THIS WAS ESTABLISHED. Not by derivation. Two hand derivations were
 * written here first and both were wrong. The sign is a parity of a few
 * crossing counts, so the candidates were enumerated instead: each of S_+ and
 * S_- may depend on any subset of {alphas below p, betas below p, alphas in
 * the determinant, betas in the determinant}, giving 16 x 16 = 256 pairs.
 * Every pair was filtered against <S^2> = S(S+1) on every CSF of thirteen
 * DRTs -- SOCI and FOCI, 2 to 5 electrons, S = 0, 1/2, 1, 3/2, 2, every Ms,
 * 1650 determinants in total, including odd-electron systems so that the two
 * totals differ in parity and Ms != 0 cases so that they differ in value.
 *
 * Exactly two pairs survived, and they differ by (n_a_total + n_b_total) on
 * BOTH operators. That is the total electron count, which neither S_+ nor
 * S_- changes, so the extra factor appears squared in S_- S_+ and cancels.
 * The two survivors are the same operator; no case can separate them and
 * none needs to. The simpler one is used.
 *
 * Note what the enumeration does and does not give. The expression contains
 * no reference to electron number or multiplicity -- it is pure fermion
 * reordering -- so a correct rule is correct for any system. But the search
 * establishes only that this rule passes the cases listed above. The
 * coverage figures are quoted so a later reader can judge that for
 * themselves rather than take the claim on trust. */
static double spin_move_sign(const uint64_t *da, const uint64_t *db, int p)
{
    int n = 0;
    for (int q = 0; q < p; q++) {
        if (bit_of(da, q)) n++;
        if (bit_of(db, q)) n++;
    }
    return (n & 1) ? -1.0 : 1.0;
}

/* One term of a determinant-basis vector. Fixed at MAXW words; see the note
 * on MAXW at the top of the file. */
typedef struct { uint64_t a[MAXW], b[MAXW]; double c; } dterm_t;

static int dterm_is(const dterm_t *t, const uint64_t *a, const uint64_t *b,
                    int nw)
{
    return words_equal(t->a, a, nw) && words_equal(t->b, b, nw);
}

/* Apply S_+ then S_- to the expansion (a, b, c) of one CSF, and return
 * <psi| S_- S_+ |psi>. Requires nwords <= MAXW; callers check. */
static double expect_sminus_splus(const uint64_t *a, const uint64_t *b,
                                  const double *c, int ndets, int norb)
{
    const int nw  = (norb + 63) / 64;
    const int cap = ndets * norb * norb + 1;
    dterm_t *t = malloc(sizeof(dterm_t) * (size_t)cap);
    if (!t) return NAN;
    int nt = 0;

    uint64_t a1[MAXW], b1[MAXW], a2[MAXW], b2[MAXW];

    for (int k = 0; k < ndets; k++) {
        if (c[k] == 0.0) continue;
        const uint64_t *ak = a + (size_t)k * nw;
        const uint64_t *bk = b + (size_t)k * nw;

        for (int p = 0; p < norb; p++) {
            /* S_+ at p: needs a beta and no alpha at p */
            if (!bit_of(bk, p)) continue;
            if (bit_of(ak, p))  continue;
            memcpy(a1, ak, sizeof(uint64_t) * (size_t)nw);
            memcpy(b1, bk, sizeof(uint64_t) * (size_t)nw);
            bit_set(a1, p);
            bit_clear(b1, p);
            const double c1 = c[k] * spin_move_sign(ak, bk, p);

            for (int q = 0; q < norb; q++) {
                /* S_- at q: needs an alpha and no beta at q */
                if (!bit_of(a1, q)) continue;
                if (bit_of(b1, q))  continue;
                memcpy(a2, a1, sizeof(uint64_t) * (size_t)nw);
                memcpy(b2, b1, sizeof(uint64_t) * (size_t)nw);
                bit_clear(a2, q);
                bit_set(b2, q);
                const double c2 = c1 * spin_move_sign(a1, b1, q);

                int found = -1;
                for (int i = 0; i < nt; i++)
                    if (dterm_is(&t[i], a2, b2, nw)) { found = i; break; }
                if (found >= 0) {
                    t[found].c += c2;
                } else if (nt < cap) {
                    memcpy(t[nt].a, a2, sizeof(uint64_t) * (size_t)nw);
                    memcpy(t[nt].b, b2, sizeof(uint64_t) * (size_t)nw);
                    t[nt].c = c2;
                    nt++;
                }
            }
        }
    }

    double s = 0.0;
    for (int k = 0; k < ndets; k++) {
        if (c[k] == 0.0) continue;
        const uint64_t *ak = a + (size_t)k * nw;
        const uint64_t *bk = b + (size_t)k * nw;
        for (int i = 0; i < nt; i++)
            if (dterm_is(&t[i], ak, bk, nw)) { s += c[k] * t[i].c; break; }
    }
    free(t);
    return s;
}

/* <S^2> = <S_- S_+> + Ms^2 + Ms, for every CSF of a DRT at every Ms. */
static void test_spin_purity(const char *label, hexpr_drt_kind_t kind,
                             int nve, int spin_x2, int ncore, int nval, int next)
{
    printf("-- S^2 %s: kind=%d nve=%d spin_x2=%d core=%d val=%d ext=%d --\n",
           label, (int)kind, nve, spin_x2, ncore, nval, next);

    hexpr_drt_t *drt = hexpr_drt_create(kind, nve, spin_x2, ncore, nval, next);
    if (!drt) { printf("  SKIP: create returned NULL (%s)\n", hexpr_last_error()); return; }

    const int norb  = hexpr_drt_norb(drt);
    const int ncsfs = hexpr_drt_ncsfs(drt);
    if ((norb + 63) / 64 > MAXW) {
        printf("  SKIP: norb=%d needs more than %d words; csf2det handles it, "
               "this test does not (see the MAXW note)\n", norb, MAXW);
        hexpr_drt_destroy(drt);
        return;
    }

    const double S    = 0.5 * (double)spin_x2;
    const double want = S * (S + 1.0);

    for (int ms = -spin_x2; ms <= spin_x2; ms += 2) {
        const double Ms = 0.5 * (double)ms;
        double worst = 0.0;
        int bad = 0;

        for (int i = 0; i < ncsfs; i++) {
            uint64_t *a, *b; double *c;
            int nd = expand(drt, i, ms, &a, &b, &c);
            if (nd < 0) { bad++; continue; }

            double sp = expect_sminus_splus(a, b, c, nd, norb);
            free(a); free(b); free(c);
            if (isnan(sp)) { bad++; continue; }

            double err = fabs(sp + Ms * Ms + Ms - want);
            if (err > worst) worst = err;
            if (err > 1e-10) bad++;
        }

        char buf[128];
        snprintf(buf, sizeof buf,
                 "ms_x2=%+d  <S^2> == %.2f for all %d CSFs (worst %.1e)",
                 ms, want, ncsfs, worst);
        ok(bad == 0, buf);
    }

    hexpr_drt_destroy(drt);
}

/* ------------------------------------------------------------------ */
/* (3) orthonormality and (4) Sz, over every Ms                        */
/* ------------------------------------------------------------------ */
static void test_drt(const char *label, hexpr_drt_kind_t kind,
                     int nve, int spin_x2, int ncore, int nval, int next)
{
    printf("-- (3)(4) %s: kind=%d nve=%d spin_x2=%d core=%d val=%d ext=%d --\n",
           label, (int)kind, nve, spin_x2, ncore, nval, next);

    hexpr_drt_t *drt = hexpr_drt_create(kind, nve, spin_x2, ncore, nval, next);
    if (!drt) { printf("  SKIP: create returned NULL (%s)\n", hexpr_last_error()); return; }

    const int norb   = hexpr_drt_norb(drt);
    const int ncsfs  = hexpr_drt_ncsfs(drt);
    const int nwords = (norb + 63) / 64;

    for (int ms = -spin_x2; ms <= spin_x2; ms += 2) {
        int sz_bad = 0, norm_bad = 0, orth_bad = 0, expand_bad = 0;
        double worst_norm = 0.0, worst_orth = 0.0;

        for (int i = 0; i < ncsfs; i++) {
            uint64_t *a, *b; double *c;
            int nd = expand(drt, i, ms, &a, &b, &c);
            if (nd < 0) { expand_bad++; continue; }

            for (int k = 0; k < nd; k++) {
                int na = popcount_words(a + (size_t)k * nwords, nwords);
                int nb = popcount_words(b + (size_t)k * nwords, nwords);
                if (na - nb != ms) sz_bad++;
            }

            double nrm = 0.0;
            for (int k = 0; k < nd; k++) nrm += c[k] * c[k];
            if (fabs(nrm - 1.0) > TOL) norm_bad++;
            if (fabs(nrm - 1.0) > worst_norm) worst_norm = fabs(nrm - 1.0);

            free(a); free(b); free(c);
        }

        for (int i = 0; i < ncsfs; i++) {
            for (int j = i + 1; j < ncsfs; j++) {
                int bad = 0;
                double s = overlap(drt, i, j, ms, &bad);
                if (bad) { expand_bad++; continue; }
                if (fabs(s) > TOL) orth_bad++;
                if (fabs(s) > worst_orth) worst_orth = fabs(s);
            }
        }

        char buf[128];
        snprintf(buf, sizeof buf,
                 "ms_x2=%+d  Sz ok, |norm-1|<=%.1e, |<I|J>|<=%.1e",
                 ms, worst_norm, worst_orth);
        ok(sz_bad == 0 && norm_bad == 0 && orth_bad == 0 && expand_bad == 0, buf);
    }

    hexpr_drt_destroy(drt);
}

/* Report the packing a DRT actually uses, so a reader can see from the log
 * that the multi-word path ran rather than taking the case names on trust. */
static void report_width(const char *label, hexpr_drt_kind_t kind,
                         int nve, int spin_x2, int ncore, int nval, int next)
{
    hexpr_drt_t *drt = hexpr_drt_create(kind, nve, spin_x2, ncore, nval, next);
    if (!drt) { printf("  %-16s SKIP (%s)\n", label, hexpr_last_error()); return; }
    const int norb = hexpr_drt_norb(drt);
    printf("  %-16s norb=%-3d nwords=%-2d ncsfs=%d\n",
           label, norb, (norb + 63) / 64, hexpr_drt_ncsfs(drt));
    hexpr_drt_destroy(drt);
}

/* ------------------------------------------------------------------ */
int main(void)
{
    if (hexpr_init() != HEXPR_OK) {
        printf("hexpr_init failed: %s\n", hexpr_last_error());
        return 1;
    }

    printf("=== csf2det ===\n");

    test_orbital_direction();
    test_singlet_pair();
    test_triplet_pair();

    test_drt("SOCI singlet", HEXPR_DRT_SOCI, 2, 0, 1, 1, 2);
    test_drt("FOCI singlet", HEXPR_DRT_FOCI, 2, 0, 1, 1, 2);
    test_drt("SOCI triplet", HEXPR_DRT_SOCI, 2, 2, 1, 1, 2);
    test_drt("FOCI triplet", HEXPR_DRT_FOCI, 2, 2, 1, 1, 2);
    test_drt("SOCI 4e S=1",  HEXPR_DRT_SOCI, 4, 2, 1, 2, 2);
    test_drt("SOCI 4e S=2",  HEXPR_DRT_SOCI, 4, 4, 1, 2, 2);
    test_drt("SOCI 3e S=1/2", HEXPR_DRT_SOCI, 3, 1, 1, 2, 2);

    /* The S^2 list mirrors the cases the sign rule was established on:
     * both kinds, even and ODD electron counts, S from 0 to 2, every Ms. */
    test_spin_purity("SOCI 2e S=0",   HEXPR_DRT_SOCI, 2, 0, 1, 1, 2);
    test_spin_purity("FOCI 2e S=0",   HEXPR_DRT_FOCI, 2, 0, 1, 1, 2);
    test_spin_purity("SOCI 2e S=1",   HEXPR_DRT_SOCI, 2, 2, 1, 1, 2);
    test_spin_purity("FOCI 2e S=1",   HEXPR_DRT_FOCI, 2, 2, 1, 1, 2);
    test_spin_purity("SOCI 4e S=0",   HEXPR_DRT_SOCI, 4, 0, 1, 2, 2);
    test_spin_purity("SOCI 4e S=1",   HEXPR_DRT_SOCI, 4, 2, 1, 2, 2);
    test_spin_purity("SOCI 4e S=2",   HEXPR_DRT_SOCI, 4, 4, 1, 2, 2);
    test_spin_purity("FOCI 4e S=1",   HEXPR_DRT_FOCI, 4, 2, 1, 2, 2);
    test_spin_purity("SOCI 3e S=1/2", HEXPR_DRT_SOCI, 3, 1, 1, 2, 2);
    test_spin_purity("SOCI 3e S=3/2", HEXPR_DRT_SOCI, 3, 3, 1, 2, 2);
    test_spin_purity("FOCI 3e S=1/2", HEXPR_DRT_FOCI, 3, 1, 1, 2, 2);
    test_spin_purity("SOCI 5e S=1/2", HEXPR_DRT_SOCI, 5, 1, 1, 2, 2);
    test_spin_purity("SOCI 5e S=3/2", HEXPR_DRT_SOCI, 5, 3, 1, 2, 2);

    /* ---------------- multi-word determinants ----------------
     * Everything above has norb <= 5, so csf2det's packing loop ran with
     * nwords = 1 throughout and the bit arithmetic never crossed a word
     * boundary. These DRTs put orbitals 64..69 in a second, partly filled
     * word. FOCI keeps the CSF counts small enough for the quadratic
     * orthonormality loop on the first two. */
    printf("-- multi-word determinants: the DRTs used below --\n");
    report_width("wide 2e S=0",   HEXPR_DRT_FOCI, 2, 0, 1, 1, 68);
    report_width("wide 2e S=1",   HEXPR_DRT_FOCI, 2, 2, 1, 1, 68);
    report_width("wide 3e S=1/2", HEXPR_DRT_FOCI, 3, 1, 1, 2, 66);
    report_width("wide 4e S=0",   HEXPR_DRT_FOCI, 4, 0, 1, 2, 66);

    test_drt("wide 2e S=0", HEXPR_DRT_FOCI, 2, 0, 1, 1, 68);
    test_drt("wide 2e S=1", HEXPR_DRT_FOCI, 2, 2, 1, 1, 68);

    test_spin_purity("wide 2e S=0",   HEXPR_DRT_FOCI, 2, 0, 1, 1, 68);
    test_spin_purity("wide 2e S=1",   HEXPR_DRT_FOCI, 2, 2, 1, 1, 68);
    test_spin_purity("wide 3e S=1/2", HEXPR_DRT_FOCI, 3, 1, 1, 2, 66);
    test_spin_purity("wide 4e S=0",   HEXPR_DRT_FOCI, 4, 0, 1, 2, 66);

    if (failures == 0) printf("ALL CSF2DET TESTS PASSED\n");
    else               printf("%d TEST(S) FAILED\n", failures);
    return failures;
}
