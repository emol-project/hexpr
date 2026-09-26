/* ====================================================================
 *  wigner_init.c -- precompute the lookup tables used by every
 *                   Wigner-coefficient routine.
 *
 *  Mathematically this builds:
 *      fac[k]      = k!                               (k = 0..56)
 *      rt[k]       = sqrt(k)                          (k = 0..56)
 *      den[k]      = sqrt(k!)                         (k = 0..56)
 *      fasq[2k+40] = sqrt(k!)                         (k = 0..56)
 *                    (all other entries left at 0)
 *      gamma[]     = Gamma(n/2 + 1) table             (n = 0..112)
 *
 *  The fasq[] layout is sparse and looks odd, but it's exactly what
 *  delta() expects: delta's index transformation (ia+ib-ic+40 in
 *  the original 1-based form, which is +39 in 0-based C) lands on
 *  the even slot that holds sqrt((ja+jb-jc)!).  Slots that don't
 *  correspond to a non-negative factorial argument stay zero, so
 *  forbidden triangles produce 0.0 automatically.
 *
 *  Behaviour is bit-for-bit identical to the original setwig_();
 *  only the surface syntax has changed.
 * ====================================================================
 */
#include <math.h>
#include <string.h>

#include "wigner.h"
#include "wigner_internal.h"

/* The single table instance.  Zero-initialised at load time so that
 * unfilled fasq[] slots stay 0 even before wigner_init() runs. */
wig_table_t wig_table;

void wigner_init(void)
{
    int k, i;

    /* Wipe fasq[] explicitly: only some slots are written below, and
     * a re-initialisation must not leave stale values behind. */
    memset(wig_table.fasq, 0, sizeof wig_table.fasq);

    /* fac[k] = k!,  rt[k] = sqrt(k),  den[k] = sqrt(k!),
     * and fasq[2k+40] = sqrt(k!) for k = 0..56.
     *
     * The sparse fasq[] layout is what delta() expects: its argument
     * transformation (ia+ib-ic+40, in 1-based form, i.e. +39 here)
     * lands on the even slot that holds sqrt((ja+jb-jc)!).
     *
     * The base offset is 40, NOT 38: the original Fortran wrote
     *     fasq( (k<<1) + 38 ) = den(k)              for k = 1..57
     * which in 0-based C corresponds to indices 40, 42, ..., 152
     * (i.e. 2m + 40 for m = 0..56).  Getting this off by two
     * silently breaks delta() for forbidden triangles. */
    wig_table.fac[0] = 1.0;
    for (k = 1; k < WIG_NFAC; ++k) {
        wig_table.fac[k] = wig_table.fac[k - 1] * (double)k;
    }
    for (k = 0; k < WIG_NFAC; ++k) {
        wig_table.den[k]            = sqrt(wig_table.fac[k]);
        wig_table.fasq[2 * k + 40]  = wig_table.den[k];
        wig_table.rt[k]             = sqrt((double)k);
    }

    /* Gamma(n/2 + 1) table, built with the same trick as the
     * original code:
     *   - integer-argument cells (even n) are obtained by repeated
     *     multiplication by k/2;
     *   - half-integer cells (odd n) are then multiplied by sqrt(pi).
     *
     * gamma[] is not used by the 6j/9j routines but is kept for
     * compatibility with downstream callers that may rely on it. */
    wig_table.gamma[0] = 1.0;
    wig_table.gamma[1] = 1.0;
    for (i = 1; i <= 112; ++i) {
        wig_table.gamma[i + 1] = wig_table.gamma[i - 1] * (double)i * 0.5;
    }
    /* sqrt(pi) = 2 * sqrt(atan(1)) = 2 * sqrt(pi/4). */
    wig_table.gamma[0] = 2.0 * sqrt(atan(1.0));
    for (i = 1; i <= 112; i += 2) {
        wig_table.gamma[i + 1] *= wig_table.gamma[0];
    }
}
