/* test_pair_one_4o4e.c -- hexpr_pair_eval 1e, s_full_4o4e */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "../include/hexpr.h"

static int failures = 0;
static int checks   = 0;

static void check(const char *label, int ok)
{
    checks++;
    printf("  %-60s %s\n", label, ok ? "PASS" : "FAIL");
    if (!ok) failures++;
}

static void check_status(const char *label, hexpr_status_t got, hexpr_status_t expected)
{
    checks++;
    int ok = (got == expected);
    printf("  %-60s %s  (got %d, expected %d)\n",
           label, ok ? "PASS" : "FAIL", (int)got, (int)expected);
    if (!ok) failures++;
}

/* parallel sort of (pqrs, ij, coef) by pqrs */
static void sort_by_pqrs(int n, int32_t *ij, int32_t *pqrs, double *coef)
{
    /* insertion sort -- small n */
    for (int i = 1; i < n; i++) {
        int32_t kp = pqrs[i], ki = ij[i]; double kc = coef[i];
        int j = i - 1;
        while (j >= 0 && pqrs[j] > kp) {
            pqrs[j+1] = pqrs[j]; ij[j+1] = ij[j]; coef[j+1] = coef[j];
            j--;
        }
        pqrs[j+1] = kp; ij[j+1] = ki; coef[j+1] = kc;
    }
}

int main(void)
{
    printf("=== pair_eval 1e test (s_full_4o4e, bra=0 ket=0) ===\n");

    hexpr_drt_t *drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 4, 0);
    check("hexpr_drt_create SOCI 4o4e succeeded", drt != NULL);
    if (!drt) { printf("FAIL: aborting\n"); return 1; }

    int ncsfs = hexpr_drt_ncsfs(drt);
    check("ncsfs == 20", ncsfs == 20);

    /* ------------------------------------------------------------------ */
    /* Test 1: pair_count for (0,0) with HEXPR_EXPR_ONE                   */
    /* ------------------------------------------------------------------ */
    printf("\n--- Test 1: hexpr_pair_count (bra=0, ket=0, ONE) ---\n");
    {
        int cnt = -1;
        hexpr_status_t st = hexpr_pair_count(drt, 0, 0, HEXPR_EXPR_ONE, &cnt);
        check_status("pair_count returns OK", st, HEXPR_OK);
        check("pair_count == 2", cnt == 2);
    }

    /* ------------------------------------------------------------------ */
    /* Test 2: pair_eval for (0,0) with HEXPR_EXPR_ONE                    */
    /* Reference (pair_one.txt): pqrs=5 coef~=2.0, pqrs=9 coef~=2.0        */
    /* ------------------------------------------------------------------ */
    printf("\n--- Test 2: hexpr_pair_eval (bra=0, ket=0, ONE) ---\n");
    {
        int32_t ij[8], pqrs[8];
        double  coef[8];
        int cnt = -1;
        hexpr_status_t st = hexpr_pair_eval(drt, 0, 0, HEXPR_EXPR_ONE,
                                             &cnt, ij, pqrs, coef);
        check_status("pair_eval returns OK", st, HEXPR_OK);
        check("pair_eval cnt == 2", cnt == 2);

        if (cnt == 2) {
            sort_by_pqrs(cnt, ij, pqrs, coef);
            /* term 0: pqrs=5, ij=triangle(0,0)=0, coef~=2.0 */
            check("term[0] ij == 0",    ij[0]   == 0);
            check("term[0] pqrs == 5",  pqrs[0] == 5);
            check("term[0] coef ~= 2.0", fabs(coef[0] - 2.0) < 1e-10);
            /* term 1: pqrs=9, ij=0, coef~=2.0 */
            check("term[1] ij == 0",    ij[1]   == 0);
            check("term[1] pqrs == 9",  pqrs[1] == 9);
            check("term[1] coef ~= 2.0", fabs(coef[1] - 2.0) < 1e-10);
        }
    }

    /* ------------------------------------------------------------------ */
    /* Test 3: invalid inputs                                              */
    /* ------------------------------------------------------------------ */
    printf("\n--- Test 3: invalid inputs ---\n");
    {
        int cnt;
        hexpr_status_t st;
        st = hexpr_pair_count(drt, -1, 0, HEXPR_EXPR_ONE, &cnt);
        check_status("pair_count bra=-1 -> RANGE", st, HEXPR_ERR_RANGE);
        st = hexpr_pair_count(drt, 0, ncsfs, HEXPR_EXPR_ONE, &cnt);
        check_status("pair_count ket=ncsfs -> RANGE", st, HEXPR_ERR_RANGE);
        st = hexpr_pair_count(NULL, 0, 0, HEXPR_EXPR_ONE, &cnt);
        check_status("pair_count NULL drt -> INVALID_ARG", st, HEXPR_ERR_INVALID_ARG);
    }

    /* ------------------------------------------------------------------ */
    /* Test 4: pair_count diagonal self-consistency                        */
    /* pair_count(i,i,ONE) must equal pair_eval count for all i           */
    /* ------------------------------------------------------------------ */
    printf("\n--- Test 4: diagonal self-consistency (all %d CSFs) ---\n", ncsfs);
    {
        int all_ok = 1;
        for (int i = 0; i < ncsfs; i++) {
            int cnt_count = -1, cnt_eval = -1;
            hexpr_pair_count(drt, i, i, HEXPR_EXPR_ONE, &cnt_count);
            hexpr_pair_eval(drt, i, i, HEXPR_EXPR_ONE, &cnt_eval,
                            NULL, NULL, NULL);
            if (cnt_count != cnt_eval) {
                printf("  FAIL: CSF %d: pair_count=%d pair_eval=%d\n",
                       i, cnt_count, cnt_eval);
                all_ok = 0; failures++; checks++;
            }
        }
        checks++;
        printf("  %-60s %s\n", "diagonal self-consistency all CSFs",
               all_ok ? "PASS" : "FAIL");
        if (!all_ok) failures++;
    }

    printf("\n");
    printf("=== pair_eval 1e: %d/%d checks passed ===\n", checks - failures, checks);
    if (failures == 0)
        printf("ALL TESTS PASSED\n");
    else
        printf("%d TEST(S) FAILED\n", failures);

    hexpr_drt_destroy(drt);
    return failures > 0 ? 1 : 0;
}
