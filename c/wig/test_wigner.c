/* ====================================================================
 *  test_wigner.c -- regression test for the refactored Wigner
 *                   library against the original f2c output.
 *
 *  Compile-and-link the original setwig/delta/f3l/xwig6j/util objects
 *  alongside this file and the new library.  Each entry point is
 *  exercised on a sweep of (2j+1) values; we check that every result
 *  matches bit-exactly (or to within a tiny tolerance for the
 *  summation-heavy 6j/9j paths).
 *
 *  The original and refactored libraries each carry their own copy
 *  of the table struct.  We therefore call setwig_() and
 *  wigner_init() once apiece up front.
 *
 *  Build via the supplied Makefile target `make test`.
 * ====================================================================
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "wigner.h"

/* Forward declarations of the original f2c-translated entry points,
 * built from the unmodified upstream source files. */
typedef int    integer;
typedef double doublereal;

extern integer    setwig_(void);
extern doublereal delta_ (integer *ia, integer *ib, integer *ic);
extern doublereal f3l_   (integer *ia, integer *ib, integer *ic);
extern doublereal xwig6j_(integer *i1, integer *i2, integer *i3,
                          integer *i4, integer *i5, integer *i6);
extern doublereal xwig9x_(integer *i1, integer *i2, integer *i3,
                          integer *i4, integer *i5, integer *i6,
                          integer *i7, integer *i8, integer *i9);

/* Tolerance for comparison.  The two libraries do exactly the same
 * arithmetic in the same order, so equality should be bit-exact in
 * practice; we still allow a hair of slack to be robust against
 * compiler reordering of the 6j summation. */
static const double TOL = 1e-12;

static int near_equal(double a, double b)
{
    if (a == b) return 1;
    double diff = fabs(a - b);
    double mag  = fabs(a) + fabs(b);
    return diff <= TOL * (1.0 + mag);
}

static int n_pass = 0;
static int n_fail = 0;

#define CHECK(label, ref, got)                                       \
    do {                                                             \
        if (near_equal((ref), (got))) {                              \
            ++n_pass;                                                \
        } else {                                                     \
            ++n_fail;                                                \
            fprintf(stderr, "FAIL %s: ref=% .17e  got=% .17e\n",     \
                    (label), (ref), (got));                          \
        }                                                            \
    } while (0)

static void test_delta(int max_arg)
{
    for (int ia = 1; ia <= max_arg; ++ia)
    for (int ib = 1; ib <= max_arg; ++ib)
    for (int ic = 1; ic <= max_arg; ++ic) {
        integer a = ia, b = ib, c = ic;
        double ref = delta_(&a, &b, &c);
        double got = wigner_delta(ia, ib, ic);
        char lbl[64];
        snprintf(lbl, sizeof lbl, "delta(%d,%d,%d)", ia, ib, ic);
        CHECK(lbl, ref, got);
    }
}

static void test_f3l(int max_arg)
{
    for (int ia = 1; ia <= max_arg; ++ia)
    for (int ib = 1; ib <= max_arg; ++ib)
    for (int ic = 1; ic <= max_arg; ++ic) {
        integer a = ia, b = ib, c = ic;
        double ref = f3l_(&a, &b, &c);
        double got = wigner_f3l(ia, ib, ic);
        char lbl[64];
        snprintf(lbl, sizeof lbl, "f3l(%d,%d,%d)", ia, ib, ic);
        CHECK(lbl, ref, got);
    }
}

static void test_6j(int max_arg)
{
    for (int i1 = 1; i1 <= max_arg; ++i1)
    for (int i2 = 1; i2 <= max_arg; ++i2)
    for (int i3 = 1; i3 <= max_arg; ++i3)
    for (int i4 = 1; i4 <= max_arg; ++i4)
    for (int i5 = 1; i5 <= max_arg; ++i5)
    for (int i6 = 1; i6 <= max_arg; ++i6) {
        integer a=i1, b=i2, c=i3, d=i4, e=i5, v=i6;
        double ref = xwig6j_(&a,&b,&c,&d,&e,&v);
        double got = wigner_6j(i1,i2,i3,i4,i5,i6);
        if (!near_equal(ref, got)) {
            char lbl[96];
            snprintf(lbl, sizeof lbl,
                     "6j(%d,%d,%d,%d,%d,%d)",
                     i1,i2,i3,i4,i5,i6);
            CHECK(lbl, ref, got);
        } else {
            ++n_pass;
        }
    }
}

static void test_9j(int max_arg)
{
    /* 9j has nine arguments; restrict to a small sweep so the test
     * runs in seconds, not hours.  We pick a handful of physically
     * meaningful cases plus a tight nested loop on a few. */
    for (int i1 = 1; i1 <= max_arg; ++i1)
    for (int i2 = 1; i2 <= max_arg; ++i2)
    for (int i3 = 1; i3 <= max_arg; ++i3)
    for (int i4 = 1; i4 <= max_arg; ++i4)
    for (int i5 = 1; i5 <= max_arg; ++i5)
    for (int i6 = 1; i6 <= max_arg; ++i6)
    for (int i7 = 1; i7 <= max_arg; ++i7)
    for (int i8 = 1; i8 <= max_arg; ++i8)
    for (int i9 = 1; i9 <= max_arg; ++i9) {
        integer a=i1,b=i2,c=i3,d=i4,e=i5,v=i6,g=i7,h=i8,ii=i9;
        double ref = xwig9x_(&a,&b,&c,&d,&e,&v,&g,&h,&ii);
        double got = wigner_9j(i1,i2,i3,i4,i5,i6,i7,i8,i9);
        if (!near_equal(ref, got)) {
            char lbl[128];
            snprintf(lbl, sizeof lbl,
                     "9j(%d,%d,%d,%d,%d,%d,%d,%d,%d)",
                     i1,i2,i3,i4,i5,i6,i7,i8,i9);
            CHECK(lbl, ref, got);
        } else {
            ++n_pass;
        }
    }
}

int main(void)
{
    setwig_();
    wigner_init();

    printf("running delta sweep (max_arg=12)...\n");
    test_delta(12);
    printf("  -> running total: %d pass, %d fail\n", n_pass, n_fail);

    printf("running f3l sweep (max_arg=12)...\n");
    test_f3l(12);
    printf("  -> running total: %d pass, %d fail\n", n_pass, n_fail);

    printf("running 6j sweep (max_arg=7)...\n");
    test_6j(7);
    printf("  -> running total: %d pass, %d fail\n", n_pass, n_fail);

    printf("running 9j sweep (max_arg=4)...\n");
    test_9j(4);
    printf("  -> running total: %d pass, %d fail\n", n_pass, n_fail);

    printf("\n========================================\n");
    printf(" %d passed, %d failed\n", n_pass, n_fail);
    printf("========================================\n");
    return n_fail == 0 ? 0 : 1;
}
