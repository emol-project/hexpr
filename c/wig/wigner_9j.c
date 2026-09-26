/* ====================================================================
 *  wigner_9j.c -- Wigner 9j symbol via the standard expansion in
 *                 6j symbols, refactored from xwig9x.c.
 *
 *      { ja jb jc }                    sum over k
 *      { jd je jv } = sum_k (-1)^(2k) (2k+1)
 *      { jg jh ji }    * { ja jb jc; jv ji k }
 *                      * { jd je jv; jb k  jh }
 *                      * { jg jh ji;  k ja jd }
 *
 *  In the (2j+1) integer convention used here, "(2k+1)" is just
 *  the integer ir, and "k" runs in steps of 1, which corresponds
 *  to ir running in steps of 2.
 *
 *  The original code accumulates the sum with a unary minus at the
 *  end and multiplies by (-1)^ks; we keep that exact phase.
 * ====================================================================
 */
#include <stdlib.h>

#include "wigner.h"
#include "wigner_internal.h"

static inline int imin3i(int a, int b, int c)
{
    int m = a < b ? a : b;
    return m < c ? m : c;
}
static inline int imax3i(int a, int b, int c)
{
    int m = a > b ? a : b;
    return m > c ? m : c;
}

double wigner_9j(int ia, int ib, int ic,
                 int id, int ie, int iv,
                 int ig, int ih, int ii)
{
    /* Bounds for the intermediate (2k+1).  These come from the three
     * triangle conditions  |ja-ji| <= k <= ja+ji  etc., translated
     * to (2j+1) integers with a +/- 1 adjustment. */
    int kl = imin3i(ia + ii, id + ih, ib + iv) - 1;
    int ks = imax3i(abs(ia - ii), abs(id - ih), abs(ib - iv)) + 1;

    if (ks > kl || (kl - ks) % 2 != 0) {
        return 0.0;
    }

    double sum = 0.0;
    for (int ir = ks; ir <= kl; ir += 2) {
        sum += (double)ir
             * wigner_6j(ia, ib, ic, iv, ii, ir)
             * wigner_6j(id, ie, iv, ib, ir, ih)
             * wigner_6j(ig, ih, ii, ir, ia, id);
    }
    /* The original applies (-1)^ks and the leading unary minus
     * here -- preserved verbatim. */
    return -sum * wig_phase(ks);
}
