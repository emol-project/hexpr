/* ====================================================================
 *  wigner.h -- Wigner 3j / 6j / 9j coefficient library
 *
 *  Refactored from f2c-translated C originating from F. Sasaki's own
 *  independent Fortran implementation. Copyright for the algorithm and
 *  its implementation remains with F. Sasaki; see LICENSE in this
 *  directory. Numerical output is bit-exact against the original.
 *  - All angular momentum arguments are passed as the integer (2j+1).
 *    For an integer j, pass (2j+1) which is odd.
 *    For a half-integer j, pass (2j+1) which is even.
 *    e.g.  j = 0   -> 1
 *          j = 1/2 -> 2
 *          j = 1   -> 3
 *          j = 3/2 -> 4
 *          ...
 *  - Internal tables support arguments up to (2j+1) <= 56
 *    (i.e. j <= 27.5).  Out-of-range calls return 0.0 and print a
 *    warning to stderr.
 *
 *  All functions are thread-unsafe in the strict sense because they
 *  share the precomputed table set up by wigner_init().  However,
 *  once wigner_init() has been called, concurrent reads of the table
 *  are safe; the library does not modify the table after init.
 * ====================================================================
 */
#ifndef WIGNER_H_
#define WIGNER_H_

#ifdef __cplusplus
extern "C" {
#endif

/* ----------------------------------------------------------------
 *  Initialization.
 *
 *  wigner_init() builds the factorial / square-root tables used by
 *  every other routine.  Call it exactly once before any other
 *  function in this header.  Calling it again simply rebuilds the
 *  tables (idempotent).
 * ---------------------------------------------------------------- */
void   wigner_init(void);

/* ----------------------------------------------------------------
 *  Triangle coefficient
 *
 *      Delta(a,b,c) = sqrt( (a+b-c)! (a-b+c)! (-a+b+c)! / (a+b+c+1)! )
 *
 *  The arguments are the (2j+1) integers, NOT j itself.
 *  The original code returns only the numerator product
 *      sqrt((a+b-c)!) * sqrt((a-b+c)!) * sqrt((-a+b+c)!)
 *  (i.e. WITHOUT the (a+b+c+1)! denominator); f3l() supplies that
 *  factor.  This behaviour is preserved verbatim because every
 *  caller in the original library relies on it.
 * ---------------------------------------------------------------- */
double wigner_delta(int two_ja_plus_1,
                    int two_jb_plus_1,
                    int two_jc_plus_1);

/* ----------------------------------------------------------------
 *  Triangle factor used in 3j evaluation
 *
 *  Returns 0 if the triangle is forbidden.
 * ---------------------------------------------------------------- */
double wigner_f3l  (int two_ja_plus_1,
                    int two_jb_plus_1,
                    int two_jc_plus_1);

/* ----------------------------------------------------------------
 *  Wigner 6j symbol with special-case handling for j == 0
 *  (i.e. argument == 1).
 *
 *      { j1 j2 j3 }
 *      { j4 j5 j6 }
 *
 *  Returns 0.0 if the symbol is forbidden by triangular conditions
 *  or if the internal table size is exceeded.
 * ---------------------------------------------------------------- */
double wigner_6j   (int two_j1_plus_1, int two_j2_plus_1, int two_j3_plus_1,
                    int two_j4_plus_1, int two_j5_plus_1, int two_j6_plus_1);

/* ----------------------------------------------------------------
 *  Wigner 9j symbol
 *
 *      { ja jb jc }
 *      { jd je jv }
 *      { jg jh ji }
 *
 *  Implemented as the standard sum over 6j symbols.
 * ---------------------------------------------------------------- */
double wigner_9j   (int two_ja_plus_1, int two_jb_plus_1, int two_jc_plus_1,
                    int two_jd_plus_1, int two_je_plus_1, int two_jv_plus_1,
                    int two_jg_plus_1, int two_jh_plus_1, int two_ji_plus_1);

#ifdef __cplusplus
}
#endif

#endif /* WIGNER_H_ */
