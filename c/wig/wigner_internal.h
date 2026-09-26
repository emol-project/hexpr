/* ====================================================================
 *  wigner_internal.h -- internal definitions for the Wigner library.
 *
 *  Not part of the public API.  Anything declared here may change
 *  between releases.
 * ====================================================================
 */
#ifndef WIGNER_INTERNAL_H_
#define WIGNER_INTERNAL_H_

/* ----------------------------------------------------------------
 *  Table sizes.
 *
 *  These match the original Fortran/f2c code 1:1 to preserve the
 *  index arithmetic that downstream routines depend on.
 *
 *      WIG_NFAC   : size of fac[] = k!  for k = 0..56
 *      WIG_NFASQ  : size of fasq[] (sparse table indexed by 2k+38)
 *      WIG_NRT    : size of rt[]  = sqrt(k) for k = 0..56
 *      WIG_NDEN   : size of den[] = sqrt(k!) for k = 0..56
 *      WIG_NGAMMA : size of gamma[] (Gamma(n/2+1) table, unused by
 *                                    6j/9j but kept for binary
 *                                    compatibility / future use)
 * ---------------------------------------------------------------- */
#define WIG_NFAC     57
#define WIG_NFASQ   154
#define WIG_NRT      57
#define WIG_NDEN     57
#define WIG_NGAMMA  114
#define WIG_NDFAC   114

/* ----------------------------------------------------------------
 *  Maximum value of the loop bound `ip` in xwig6j.  Beyond this the
 *  factorial table overflows.  Kept as a named constant so the
 *  warning message stays in sync with reality.
 * ---------------------------------------------------------------- */
#define WIG_MAX_IP   56

/* ----------------------------------------------------------------
 *  Precomputed tables.
 *
 *  Indexing convention follows the original code: index 0 here
 *  corresponds to Fortran index 1.  In particular:
 *      fac[k]  = k!                       for k = 0..56
 *      rt[k]   = sqrt(k)                  for k = 0..56
 *      den[k]  = sqrt(k!)                 for k = 0..56
 *      fasq[i] is non-zero only for even i in [40, 152]:
 *              fasq[2k+40] = sqrt(k!)     for k = 0..56
 *              The base offset 40 (not 38) is what makes delta's
 *              "ia+ib-ic+39" index hit a valid slot when the
 *              triangle is allowed.
 *
 *  The struct is exposed via a single static instance defined in
 *  wigner_init.c.  Callers should treat it as read-only after
 *  wigner_init() has run.
 * ---------------------------------------------------------------- */
typedef struct {
    double rt   [WIG_NRT];
    double fasq [WIG_NFASQ];
    double den  [WIG_NDEN];
    double fac  [WIG_NFAC];
    double gamma[WIG_NGAMMA];
    double dfac [WIG_NDFAC];
} wig_table_t;

extern wig_table_t wig_table;

/* ----------------------------------------------------------------
 *  Internal helpers (not part of the public API but shared between
 *  translation units).
 *
 *  wig_phase(n) returns (-1)^n for any (possibly negative) integer n.
 *  This replaces the original pow_ii() / c_n1 trick from f2c.
 * ---------------------------------------------------------------- */
static inline double wig_phase(int n)
{
    return (n & 1) ? -1.0 : 1.0;
}

#endif /* WIGNER_INTERNAL_H_ */
