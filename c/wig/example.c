/* ====================================================================
 *  example.c -- minimal demonstration of how to call the Wigner
 *               library from a C program.
 *
 *  Build:
 *      cc -Iinclude example.c libwigner.a -lm -o example
 *  or:
 *      cc -Iinclude example.c -L. -lwigner -lm -o example
 *
 *  Run:
 *      ./example
 * ====================================================================
 */
#include <stdio.h>
#include "wigner.h"

int main(void)
{
    /* Build the precomputed factorial / sqrt-factorial tables.
     * Must be called once before any other library function.       */
    wigner_init();

    /* Arguments are the integer (2j+1).  For example j = 1 -> 3.
     *
     *     { 1 1 1 }
     *     { 1 1 1 } = 1/6
     */
    double v6 = wigner_6j(3, 3, 3,
                          3, 3, 3);
    printf("6j {1 1 1; 1 1 1}        = %.10f   (expected  0.1666666667)\n", v6);

    /*     { 1 1 0 }
     *     { 1 1 1 } = -1/3
     */
    v6 = wigner_6j(3, 3, 1,
                   3, 3, 3);
    printf("6j {1 1 0; 1 1 1}        = %.10f   (expected -0.3333333333)\n", v6);

    /*     { 1 1 0 }
     *     { 1 1 0 }
     *     { 0 0 0 } = 1/3
     */
    double v9 = wigner_9j(3, 3, 1,
                          3, 3, 1,
                          1, 1, 1);
    printf("9j {1 1 0; 1 1 0; 0 0 0} = %.10f   (expected  0.3333333333)\n", v9);

    /* Triangle factor demos. */
    printf("delta(3,3,3)             = %.10f\n", wigner_delta(3, 3, 3));
    printf("f3l(3,3,3)               = %.10f\n", wigner_f3l  (3, 3, 3));

    return 0;
}
