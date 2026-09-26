/* ====================================================================
 *  wigner_6j.c -- Wigner 6j symbol (Racah W-coefficient form),
 *                 refactored from xwig6j.c.
 *
 *  Numerical content is preserved verbatim from the original f2c
 *  output.  The reorganisation is purely surface-level:
 *
 *      - Goto-driven dispatch on "which argument is j=0" replaced by
 *        an explicit switch (one labelled case per zero argument).
 *      - The Fortran 1-based "fac[k - 1]" pattern is reduced to a
 *        plain 0-based read by adjusting the index at definition.
 *      - The pow_ii(&c_n1, &n) trick is replaced with a small
 *        inline helper wig_phase().
 *      - max-of-N / min-of-N implemented with simple loops, so the
 *        ill-conditioned global max/min functions go away.
 *
 *  Phase conventions, special-case formulae, summation bounds and
 *  the final normalisation are byte-for-byte identical to the
 *  original.
 * ====================================================================
 */
#include <stdio.h>
#include <stdlib.h>

#include "wigner.h"
#include "wigner_internal.h"

/* Tiny local helpers -- file-scoped to avoid polluting the global
 * namespace as the original max/min did. */
static inline int imax2(int a, int b) { return a > b ? a : b; }
static inline int imin2(int a, int b) { return a < b ? a : b; }
static inline int imax4(int a, int b, int c, int d)
{
    return imax2(imax2(a, b), imax2(c, d));
}
static inline int imin3(int a, int b, int c)
{
    return imin2(imin2(a, b), c);
}

/* --------------------------------------------------------------
 *  Special-case helpers for "one argument is j=0".
 *
 *  Each of these checks the corresponding consistency constraint
 *  on the remaining arguments and, if satisfied, returns the
 *  closed-form value
 *
 *      (-1)^(s/2 + 1) / sqrt((2x+1)(2y+1))
 *
 *  where (s, x, y) depend on which argument is zero.
 *
 *  The phase exponent (s/2 + 1) is computed in integer arithmetic;
 *  s is even because the surviving triangle has integer total.
 * -------------------------------------------------------------- */

/* j1 == 0  (i1 == 1):  requires  j2 == j3  and  j5 == j6. */
static double sixj_j1_zero(int ib, int ic, int id, int ie, int iv)
{
    if (ib != ic || ie != iv)        return 0.0;
    if (wigner_delta(id, ie, ic) == 0.0) return 0.0;
    int s = (id + ie + ic) / 2 + 1;
    return wig_phase(s) / (wig_table.rt[ib] * wig_table.rt[ie]);
}
/* j2 == 0  (i2 == 1):  requires  j1 == j3  and  j4 == j6. */
static double sixj_j2_zero(int ia, int ic, int id, int ie, int iv)
{
    if (ia != ic || id != iv)        return 0.0;
    if (wigner_delta(id, ie, ic) == 0.0) return 0.0;
    int s = (id + ie + ic) / 2 + 1;
    return wig_phase(s) / (wig_table.rt[ia] * wig_table.rt[id]);
}
/* j3 == 0  (i3 == 1):  requires  j1 == j2  and  j4 == j5. */
static double sixj_j3_zero(int ia, int ib, int id, int ie, int iv)
{
    if (ia != ib || id != ie)        return 0.0;
    if (wigner_delta(ia, ie, iv) == 0.0) return 0.0;
    int s = (ia + ie + iv) / 2 + 1;
    return wig_phase(s) / (wig_table.rt[ia] * wig_table.rt[id]);
}
/* j4 == 0  (i4 == 1):  requires  j6 == j2  and  j3 == j5. */
static double sixj_j4_zero(int ia, int ib, int ic, int ie, int iv)
{
    if (iv != ib || ic != ie)        return 0.0;
    if (wigner_delta(ia, ib, ic) == 0.0) return 0.0;
    int s = (ia + ib + ic) / 2 + 1;
    return wig_phase(s) / (wig_table.rt[ib] * wig_table.rt[ic]);
}
/* j5 == 0  (i5 == 1):  requires  j6 == j1  and  j3 == j4. */
static double sixj_j5_zero(int ia, int ib, int ic, int id, int iv)
{
    if (iv != ia || ic != id)        return 0.0;
    if (wigner_delta(ia, ib, ic) == 0.0) return 0.0;
    int s = (ia + ib + ic) / 2 + 1;
    return wig_phase(s) / (wig_table.rt[ia] * wig_table.rt[ic]);
}
/* j6 == 0  (i6 == 1):  requires  j5 == j1  and  j2 == j4. */
static double sixj_j6_zero(int ia, int ib, int ic, int id, int ie)
{
    if (ie != ia || ib != id)        return 0.0;
    if (wigner_delta(ia, ib, ic) == 0.0) return 0.0;
    int s = (ia + ib + ic) / 2 + 1;
    return wig_phase(s) / (wig_table.rt[ia] * wig_table.rt[ib]);
}

/* --------------------------------------------------------------
 *  Generic 6j evaluator (Racah's sum form).
 *
 *  All triangle deltas are required to be positive; otherwise zero.
 *  The summation index ir runs from
 *
 *      iq = max(ik, il, im, in)
 *
 *  to
 *
 *      ip = min(kk, kl, km)
 *
 *  with each integer factorial argument guaranteed in [0, ip]
 *  by construction.  If ip exceeds the precomputed factorial
 *  table size (WIG_MAX_IP) we abort to zero with a warning, as
 *  the original did.
 * -------------------------------------------------------------- */
static double sixj_general(int ia, int ib, int ic, int id, int ie, int iv)
{
    double a = wigner_delta(ia, ib, ic) *
               wigner_delta(id, ie, ic) *
               wigner_delta(id, ib, iv) *
               wigner_delta(ia, ie, iv);
    if (a <= 0.0) {
        return 0.0;
    }

    int ik = (ia + ib + ic) / 2;
    int il = (id + ie + ic) / 2;
    int im = (id + ib + iv) / 2;
    int in = (ia + ie + iv) / 2;

    int kk = (ia + ib + id + ie - 2) / 2;
    int kl = (ia + ic + id + iv - 2) / 2;
    int km = (ib + ic + ie + iv - 2) / 2;

    int iq = imax4(ik, il, im, in);
    int ip = imin3(kk, kl, km);

    if (ip > WIG_MAX_IP) {
        fprintf(stderr,
                "/// ***** wigner_6j RANGE EXCEEDS %d  ip=%5d  PAR HAM = "
                "%3d%3d%3d%3d%3d%3d ***** ///\n",
                WIG_MAX_IP, ip, ia, ib, ic, id, ie, iv);
        return 0.0;
    }

    double g = 0.0;
    for (int ir = iq; ir <= ip; ++ir) {
        /* Original arithmetic kept identical: jk = ir-ik+1, then
         * read fac[jk-1] = fac[ir-ik].  The "+1" / "-1" pair
         * cancels; we just use ir-ik etc. directly. */
        int jk = ir - ik;
        int jl = ir - il;
        int jm = ir - im;
        int jn = ir - in;
        int lk = kk - ir;
        int ll = kl - ir;
        int lm = km - ir;

        double denom = wig_table.fac[jk] * wig_table.fac[jl] *
                       wig_table.fac[jm] * wig_table.fac[jn] *
                       wig_table.fac[lk] * wig_table.fac[ll] *
                       wig_table.fac[lm];
        g -= wig_phase(ir) * wig_table.fac[ir] / denom;
    }

    return g * a /
           (wig_table.den[ik] * wig_table.den[il] *
            wig_table.den[im] * wig_table.den[in]);
}

/* --------------------------------------------------------------
 *  Public entry point.
 *
 *  Dispatch order matters: the original code checks i1 first, then
 *  i2, etc., and falls through to the general routine only when
 *  none of the arguments is 1.  We keep the same precedence even
 *  though the special-case formulae are mutually consistent, just
 *  to make the refactor a true behavioural no-op.
 * -------------------------------------------------------------- */
double wigner_6j(int i1, int i2, int i3, int i4, int i5, int i6)
{
    if (i1 == 1) return sixj_j1_zero(i2, i3, i4, i5, i6);
    if (i2 == 1) return sixj_j2_zero(i1, i3, i4, i5, i6);
    if (i3 == 1) return sixj_j3_zero(i1, i2, i4, i5, i6);
    if (i4 == 1) return sixj_j4_zero(i1, i2, i3, i5, i6);
    if (i5 == 1) return sixj_j5_zero(i1, i2, i3, i4, i6);
    if (i6 == 1) return sixj_j6_zero(i1, i2, i3, i4, i5);
    return sixj_general(i1, i2, i3, i4, i5, i6);
}
