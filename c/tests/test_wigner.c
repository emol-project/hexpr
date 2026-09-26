/* Smoke test: numerical spot checks on libwigner (wigner_init, wigner_6j,
 * wigner_9j).  It links ../libwigner.a directly.  It used to reach these
 * through libhexpr.so, which re-exported them; that is not a property
 * worth testing -- libhexpr.so cannot be linked at all without the wigner
 * objects, because calc_TensorOp_average calls them.  The deeper check
 * against the original f2c implementation lives in wig/test_wigner.c. */

#include <stdio.h>
#include <math.h>
#include "wigner.h"

static int check(const char *label, double got, double expected, double tol)
{
    int ok = fabs(got - expected) <= tol;
    printf("  %-50s %s  (got %.15f, expected %.15f)\n",
           label, ok ? "PASS" : "FAIL", got, expected);
    return ok ? 0 : 1;
}

int main(void)
{
    int failures = 0;
    printf("=== Wigner smoke test ===\n");

    wigner_init();

    /* All j=0 (all args = 2j+1 = 1): {0 0 0; 0 0 0; 0 0 0} = 1 */
    failures += check("wigner_9j(all-1) = 1.0",
                      wigner_9j(1,1,1, 1,1,1, 1,1,1), 1.0, 1e-14);

    /* All j=0 for 6j: {0 0 0; 0 0 0} = 1 */
    failures += check("wigner_6j(all-1) = 1.0",
                      wigner_6j(1,1,1, 1,1,1), 1.0, 1e-14);

    /* UNITY segment on a singlet--singlet path:
     * hmlkva(EMPTY, EMPTY, UNITY, psi=1, phi=1, theta=1, theta_prev=1)
     * = (-1)^0 * matrix_el(1,1,1)=1 * sqrt(1*1*1)=1 * wigner_9j(1,1,1,1,1,1,1,1,1)
     * = 1.0 (consistency check, same call as above) */
    failures += check("wigner_9j for UNITY/singlet/singlet = 1.0",
                      wigner_9j(1,1,1, 1,1,1, 1,1,1), 1.0, 1e-14);

    printf("%s (%d failure(s))\n",
           failures == 0 ? "PASS" : "FAIL", failures);
    return failures ? 1 : 0;
}
