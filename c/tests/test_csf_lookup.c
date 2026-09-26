/* test_csf_lookup.c -- hexpr_csf_index, hexpr_csf_steps, hexpr_step_parse */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

int main(void)
{
    printf("=== CSF lookup test (s_full_4o4e) ===\n");

    hexpr_drt_t *drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 4, 0);
    check("hexpr_drt_create SOCI 4o4e succeeded", drt != NULL);
    if (!drt) {
        printf("FAIL: DRT creation failed; aborting\n");
        return 1;
    }

    int norb  = hexpr_drt_norb(drt);
    int ncsfs = hexpr_drt_ncsfs(drt);
    check("norb == 4",   norb  == 4);
    check("ncsfs == 20", ncsfs == 20);

    const unsigned char *csfs = hexpr_drt_csfs(drt);
    check("hexpr_drt_csfs not NULL", csfs != NULL);
    if (!csfs) {
        printf("FAIL: csfs NULL; aborting\n");
        hexpr_drt_destroy(drt);
        return 1;
    }

    /* ------------------------------------------------------------------ */
    /* Test 1: hexpr_csf_index round-trip on all 20 CSFs                  */
    /* ------------------------------------------------------------------ */
    printf("\n--- Test 1: csf_index round-trip (all %d CSFs) ---\n", ncsfs);
    {
        int all_ok = 1;
        for (int i = 0; i < ncsfs; i++) {
            const unsigned char *step = csfs + i * norb;
            int idx = -1;
            hexpr_status_t st = hexpr_csf_index(drt, step, &idx);
            if (st != HEXPR_OK || idx != i) {
                printf("  FAIL: csf_index for CSF %d: status=%d idx=%d\n", i, (int)st, idx);
                all_ok = 0;
                failures++;
                checks++;
            }
        }
        checks++;
        printf("  %-60s %s\n",
               "csf_index round-trip all CSFs", all_ok ? "PASS" : "FAIL");
        if (!all_ok) failures++;
    }

    /* ------------------------------------------------------------------ */
    /* Test 2: hexpr_csf_steps round-trip on all 20 CSFs                  */
    /* ------------------------------------------------------------------ */
    printf("\n--- Test 2: csf_steps round-trip (all %d CSFs) ---\n", ncsfs);
    {
        unsigned char *buf = malloc((size_t)norb);
        int all_ok = 1;
        for (int i = 0; i < ncsfs; i++) {
            const unsigned char *expected = csfs + i * norb;
            hexpr_status_t st = hexpr_csf_steps(drt, i, buf);
            if (st != HEXPR_OK || memcmp(buf, expected, (size_t)norb) != 0) {
                printf("  FAIL: csf_steps for CSF %d: status=%d\n", i, (int)st);
                all_ok = 0;
                failures++;
                checks++;
            }
        }
        free(buf);
        checks++;
        printf("  %-60s %s\n",
               "csf_steps round-trip all CSFs", all_ok ? "PASS" : "FAIL");
        if (!all_ok) failures++;
    }

    /* ------------------------------------------------------------------ */
    /* Test 3: hexpr_step_parse                                            */
    /* ------------------------------------------------------------------ */
    printf("\n--- Test 3: hexpr_step_parse ---\n");
    {
        unsigned char buf[4];
        hexpr_status_t st;

        st = hexpr_step_parse("eeff", 4, buf);
        check_status("step_parse 'eeff' returns OK", st, HEXPR_OK);
        check("step_parse 'eeff' [0]==E", buf[0] == HEXPR_STEP_E);
        check("step_parse 'eeff' [1]==E", buf[1] == HEXPR_STEP_E);
        check("step_parse 'eeff' [2]==F", buf[2] == HEXPR_STEP_F);
        check("step_parse 'eeff' [3]==F", buf[3] == HEXPR_STEP_F);

        st = hexpr_step_parse("UUDD", 4, buf);
        check_status("step_parse 'UUDD' returns OK", st, HEXPR_OK);
        check("step_parse 'UUDD' [0]==U", buf[0] == HEXPR_STEP_U);
        check("step_parse 'UUDD' [1]==U", buf[1] == HEXPR_STEP_U);
        check("step_parse 'UUDD' [2]==D", buf[2] == HEXPR_STEP_D);
        check("step_parse 'UUDD' [3]==D", buf[3] == HEXPR_STEP_D);

        st = hexpr_step_parse("UdUd", 4, buf);
        check_status("step_parse 'UdUd' returns OK", st, HEXPR_OK);
        check("step_parse 'UdUd' [0]==U", buf[0] == HEXPR_STEP_U);
        check("step_parse 'UdUd' [1]==D", buf[1] == HEXPR_STEP_D);
        check("step_parse 'UdUd' [2]==U", buf[2] == HEXPR_STEP_U);
        check("step_parse 'UdUd' [3]==D", buf[3] == HEXPR_STEP_D);
    }

    /* ------------------------------------------------------------------ */
    /* Test 4: invalid inputs                                              */
    /* ------------------------------------------------------------------ */
    printf("\n--- Test 4: invalid inputs ---\n");
    {
        unsigned char step_bad[4] = {0, 0, 5, 0};  /* 5 is invalid */
        int idx = -1;
        hexpr_status_t st;

        st = hexpr_csf_index(drt, step_bad, &idx);
        check_status("csf_index with step[2]=5 returns INVALID_ARG",
                     st, HEXPR_ERR_INVALID_ARG);

        unsigned char buf[4];
        st = hexpr_csf_steps(drt, -1, buf);
        check_status("csf_steps index=-1 returns RANGE", st, HEXPR_ERR_RANGE);

        st = hexpr_csf_steps(drt, ncsfs, buf);
        check_status("csf_steps index=ncsfs returns RANGE", st, HEXPR_ERR_RANGE);

        st = hexpr_step_parse("eexx", 4, buf);
        check_status("step_parse 'eexx' returns INVALID_ARG",
                     st, HEXPR_ERR_INVALID_ARG);

        st = hexpr_step_parse(NULL, 4, buf);
        check_status("step_parse NULL str returns INVALID_ARG",
                     st, HEXPR_ERR_INVALID_ARG);

        /* step sequence valid syntax but invalid DRT path (all-'d' from root) */
        unsigned char step_nodrt[4] = {HEXPR_STEP_D, 0, 0, 0};
        st = hexpr_csf_index(drt, step_nodrt, &idx);
        check_status("csf_index with invalid DRT path returns INVALID_ARG",
                     st, HEXPR_ERR_INVALID_ARG);
    }

    /* ------------------------------------------------------------------ */
    /* Summary                                                             */
    /* ------------------------------------------------------------------ */
    printf("\n");
    printf("=== CSF lookup: %d/%d checks passed ===\n", checks - failures, checks);
    if (failures == 0)
        printf("ALL TESTS PASSED\n");
    else
        printf("%d TEST(S) FAILED\n", failures);

    hexpr_drt_destroy(drt);
    return failures > 0 ? 1 : 0;
}
