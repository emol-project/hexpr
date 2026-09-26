/* ====================================================================
 *  wigner_3j.c -- triangle coefficient (delta) and 3j building block
 *                 (f3l), refactored from delta.c / f3l.c.
 *
 *  The original f2c-translated code used 1-based indexing internally
 *  (e.g. fasq[i1 - 1]).  Here those offsets are folded into the
 *  index expression so that the array is read in plain C 0-based
 *  fashion.  The arithmetic is otherwise identical.
 *
 *  All public functions accept the standard library convention:
 *  arguments are the integer (2j+1).
 * ====================================================================
 */
#include "wigner.h"
#include "wigner_internal.h"

/* --------------------------------------------------------------
 *  wigner_delta
 *
 *  Returns the triangle factor without the (a+b+c+1)! denominator
 *  (that piece is supplied by wigner_f3l below).  Specifically:
 *
 *      delta = sqrt((ja+jb-jc)!) * sqrt((ja-jb+jc)!) *
 *              sqrt((-ja+jb+jc)!)
 *
 *  where ja = (ia-1)/2 etc.  When any of those factorial arguments
 *  is negative, the corresponding fasq[] slot is zero, and the
 *  product is therefore zero.
 *
 *  The index expressions ia+ib-ic+39 etc. are equivalent to:
 *      i = 2*(ja+jb-jc) + 38
 *  which is exactly the slot wigner_init() filled with
 *  sqrt((ja+jb-jc)!).
 * -------------------------------------------------------------- */
double wigner_delta(int ia, int ib, int ic)
{
    int i1 = ia + ib - ic + 39;
    int i2 = ia - ib + ic + 39;
    int i3 = -ia + ib + ic + 39;

    /* Defensive bounds check: out-of-range indices return 0, matching
     * the implicit behaviour of the original which relied on fasq[]
     * being zero outside the populated region. */
    if (i1 < 0 || i1 >= WIG_NFASQ ||
        i2 < 0 || i2 >= WIG_NFASQ ||
        i3 < 0 || i3 >= WIG_NFASQ) {
        return 0.0;
    }
    return wig_table.fasq[i1] * wig_table.fasq[i2] * wig_table.fasq[i3];
}

/* --------------------------------------------------------------
 *  wigner_f3l
 *
 *  Triangle factor used in 3j evaluation.  Returns 0 immediately
 *  when delta is zero (forbidden triangle).
 *
 *  Index map (in original 1-based form):
 *      it = (ia+ib+ic+1) / 4
 *      iu =  ia+ib+ic   + 40
 *      i1 = (ia+ib-ic+3) / 4
 *      i2 = (ib+ic-ia+3) / 4
 *      i3 = (ic+ia-ib+3) / 4
 *
 *  Translated to 0-based by subtracting 1 wherever the original
 *  wrote arr[idx - 1].
 *
 *  Note: integer division here is *truncated*, exactly as in the
 *  original Fortran/f2c.  We keep that behaviour because changing
 *  it would alter results.
 * -------------------------------------------------------------- */
double wigner_f3l(int ia, int ib, int ic)
{
    double d = wigner_delta(ia, ib, ic);
    if (d == 0.0) {
        return 0.0;
    }

    int it = (ia + ib + ic + 1) / 4;
    int iu =  ia + ib + ic + 40;
    int i1 = (ia + ib - ic + 3) / 4;
    int i2 = (ib + ic - ia + 3) / 4;
    int i3 = (ic + ia - ib + 3) / 4;

    /* Convert to 0-based offsets. */
    int it0 = it - 1;
    int iu0 = iu - 1;
    int i10 = i1 - 1;
    int i20 = i2 - 1;
    int i30 = i3 - 1;

    if (it0 < 0 || it0 >= WIG_NFAC  ||
        iu0 < 0 || iu0 >= WIG_NFASQ ||
        i10 < 0 || i10 >= WIG_NFAC  ||
        i20 < 0 || i20 >= WIG_NFAC  ||
        i30 < 0 || i30 >= WIG_NFAC) {
        return 0.0;
    }

    double denom = wig_table.fasq[iu0] *
                   wig_table.fac [i10] *
                   wig_table.fac [i20] *
                   wig_table.fac [i30];
    if (denom == 0.0) {
        return 0.0;
    }
    return d * wig_table.fac[it0] / denom;
}
