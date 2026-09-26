/* test_max_terms.c -- max_terms helpers */
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "../include/hexpr.h"

static int failures = 0;
static int checks   = 0;

static void check(const char *label, int ok)
{
    checks++;
    if (!ok) { failures++; printf("  %-55s FAIL\n", label); }
    else     {              printf("  %-55s PASS\n", label); }
}

int main(void)
{
    hexpr_init();
    hexpr_drt_t *drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 4, 0);
    int norb = hexpr_drt_norb(drt);
    int n_oel = norb * (norb + 1) / 2;  /* 10 */
    int per_pair_full_expected = n_oel + n_oel * (n_oel + 1) / 2;  /* 65 */
    int per_pair_one_expected  = n_oel;  /* 10 */
    int max_val;

    printf("=== max_terms test (s_full_4o4e) ===\n");
    printf("  norb=%d, n_oel=%d\n", norb, n_oel);

    /* --- Test 1: hexpr_pair_max_terms --- */
    printf("\n--- Test 1: hexpr_pair_max_terms ---\n");
    check("pair FULL returns OK",
        hexpr_pair_max_terms(drt, HEXPR_EXPR_FULL, &max_val) == HEXPR_OK);
    check("pair FULL == 65", max_val == per_pair_full_expected);
    check("pair ONE returns OK",
        hexpr_pair_max_terms(drt, HEXPR_EXPR_ONE, &max_val) == HEXPR_OK);
    check("pair ONE == 10", max_val == per_pair_one_expected);

    /* Compare actual max to upper bound across all 400 pairs */
    int ncsfs = hexpr_drt_ncsfs(drt);
    int actual_max_full = 0;
    int actual_max_one  = 0;
    for (int i = 0; i < ncsfs; i++) {
        for (int j = 0; j < ncsfs; j++) {
            int n;
            hexpr_pair_count(drt, i, j, HEXPR_EXPR_FULL, &n);
            if (n > actual_max_full) actual_max_full = n;
            hexpr_pair_count(drt, i, j, HEXPR_EXPR_ONE, &n);
            if (n > actual_max_one) actual_max_one = n;
        }
    }
    printf("  measured max FULL: %d (bound %d)\n", actual_max_full, per_pair_full_expected);
    printf("  measured max ONE:  %d (bound %d)\n",  actual_max_one, per_pair_one_expected);
    check("actual max FULL <= bound", actual_max_full <= per_pair_full_expected);
    check("actual max ONE <= bound",  actual_max_one  <= per_pair_one_expected);

    /* --- Test 2: hexpr_block_max_terms --- */
    printf("\n--- Test 2: hexpr_block_max_terms ---\n");
    check("block 1x1 FULL == 65",
        hexpr_block_max_terms(drt, 1, 1, HEXPR_EXPR_FULL, &max_val) == HEXPR_OK
        && max_val == per_pair_full_expected);
    check("block 2x3 FULL == 6*65 = 390",
        hexpr_block_max_terms(drt, 2, 3, HEXPR_EXPR_FULL, &max_val) == HEXPR_OK
        && max_val == 2 * 3 * per_pair_full_expected);
    check("block 0x5 returns 0",
        hexpr_block_max_terms(drt, 0, 5, HEXPR_EXPR_FULL, &max_val) == HEXPR_OK
        && max_val == 0);
    check("block 5x0 returns 0",
        hexpr_block_max_terms(drt, 5, 0, HEXPR_EXPR_FULL, &max_val) == HEXPR_OK
        && max_val == 0);

    /* --- Test 3: hexpr_subspace_max_terms --- */
    printf("\n--- Test 3: hexpr_subspace_max_terms ---\n");
    check("subspace n=20 FULL == 20*20*65",
        hexpr_subspace_max_terms(drt, 20, HEXPR_EXPR_FULL, &max_val) == HEXPR_OK
        && max_val == 20 * 20 * per_pair_full_expected);
    check("subspace n=0 returns 0",
        hexpr_subspace_max_terms(drt, 0, HEXPR_EXPR_FULL, &max_val) == HEXPR_OK
        && max_val == 0);

    /* --- Test 4: invalid inputs --- */
    printf("\n--- Test 4: invalid inputs ---\n");
    check("pair NULL drt -> INVALID_ARG",
        hexpr_pair_max_terms(NULL, HEXPR_EXPR_FULL, &max_val) == HEXPR_ERR_INVALID_ARG);
    check("pair NULL out_max -> INVALID_ARG",
        hexpr_pair_max_terms(drt, HEXPR_EXPR_FULL, NULL) == HEXPR_ERR_INVALID_ARG);
    check("block n_bra=-1 -> INVALID_ARG",
        hexpr_block_max_terms(drt, -1, 5, HEXPR_EXPR_FULL, &max_val) == HEXPR_ERR_INVALID_ARG);
    check("block n_ket=-1 -> INVALID_ARG",
        hexpr_block_max_terms(drt, 5, -1, HEXPR_EXPR_FULL, &max_val) == HEXPR_ERR_INVALID_ARG);
    check("subspace n=-1 -> INVALID_ARG",
        hexpr_subspace_max_terms(drt, -1, HEXPR_EXPR_FULL, &max_val) == HEXPR_ERR_INVALID_ARG);

    /* --- Test 5: integer overflow -> ERR_RANGE --- */
    printf("\n--- Test 5: integer overflow ---\n");
    /* For 4o4e (per_pair=65), block n_bra * n_ket * 65 > INT_MAX
     * when n_bra * n_ket > 33028787 (~5747 squared). */
    check("block huge n -> ERR_RANGE",
        hexpr_block_max_terms(drt, 100000, 100000, HEXPR_EXPR_FULL, &max_val) == HEXPR_ERR_RANGE);
    check("subspace huge n -> ERR_RANGE",
        hexpr_subspace_max_terms(drt, 100000, HEXPR_EXPR_FULL, &max_val) == HEXPR_ERR_RANGE);

    /* --- Summary --- */
    printf("\n=== max_terms: %d/%d checks passed ===\n", checks - failures, checks);
    if (failures > 0) {
        printf("FAIL: %d failure(s)\n", failures);
    } else {
        printf("ALL TESTS PASSED\n");
    }

    hexpr_drt_destroy(drt);
    return failures > 0 ? 1 : 0;
}
