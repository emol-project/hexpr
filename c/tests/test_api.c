/* test_api.c -- exercise every public function on s_full_4o4e */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include "../include/hexpr.h"
#include "util/compare.h"

#define REF_EXPR "../../reference/s_full_4o4e.expr.txt"
#define TMP_DIR  "tmp/test_api"

static int failures = 0;
static int checks   = 0;

static void check_int(const char *label, int got, int expected)
{
    checks++;
    int ok = (got == expected);
    printf("  %-55s %s  (got %d, expected %d)\n",
           label, ok ? "PASS" : "FAIL", got, expected);
    if (!ok) failures++;
}

static void check_ok(const char *label, hexpr_status_t st)
{
    checks++;
    int ok = (st == HEXPR_OK);
    printf("  %-55s %s  (status %d)\n", label, ok ? "PASS" : "FAIL", (int)st);
    if (!ok) failures++;
}

static void check_notnull(const char *label, const void *p)
{
    checks++;
    int ok = (p != NULL);
    printf("  %-55s %s\n", label, ok ? "PASS" : "FAIL");
    if (!ok) failures++;
}

int main(void)
{
    printf("=== Public API test (s_full_4o4e) ===\n");

    /* --- hexpr_init --- */
    check_ok("hexpr_init() returns HEXPR_OK", hexpr_init());
    check_ok("hexpr_init() is idempotent",    hexpr_init());

    /* --- version --- */
    check_int("hexpr_version_major() == 1", hexpr_version_major(), 1);
    check_int("hexpr_version_minor() == 0", hexpr_version_minor(), 0);
    check_int("hexpr_version_patch() == 0", hexpr_version_patch(), 0);
    {
        const char *vs = hexpr_version_string();
        int ok = vs && strcmp(vs, "1.0.0") == 0;
        checks++;
        printf("  %-55s %s  (got '%s')\n",
               "hexpr_version_string() == '1.0.0'", ok ? "PASS" : "FAIL", vs ? vs : "(null)");
        if (!ok) failures++;
    }

    /* --- hexpr_last_error() after success returns empty --- */
    {
        const char *err = hexpr_last_error();
        checks++;
        printf("  %-55s %s\n",
               "hexpr_last_error() not NULL", err != NULL ? "PASS" : "FAIL");
        if (!err) failures++;
    }

    /* --- DRT creation --- */
    hexpr_drt_t *drt = hexpr_drt_create(HEXPR_DRT_SOCI, 4, 0, 0, 4, 0);
    check_notnull("hexpr_drt_create SOCI 4o4e", drt);
    if (!drt) { printf("Cannot continue without DRT\n"); return 1; }

    /* --- DRT queries --- */
    check_int("hexpr_drt_kind == HEXPR_DRT_SOCI", (int)hexpr_drt_kind(drt), (int)HEXPR_DRT_SOCI);
    check_int("hexpr_drt_norb == 4",   hexpr_drt_norb(drt),   4);
    check_int("hexpr_drt_ncsfs == 20", hexpr_drt_ncsfs(drt),  20);  /* known from reference */

    /* --- hexpr_drt_csfs: check step codes match reference CSF list ---
     * Reference (s_full_4o4e.expr.txt) first CSF: ["e","e","f","f"] -> 0,0,3,3
     * Second CSF: ["e","u","d","f"] -> 0,1,2,3 */
    const unsigned char *csfs = hexpr_drt_csfs(drt);
    check_notnull("hexpr_drt_csfs not NULL", csfs);
    if (csfs) {
        int norb = hexpr_drt_norb(drt);
        /* CSF 0: e,e,f,f */
        int ok0 = (csfs[0*norb+0] == 0 && csfs[0*norb+1] == 0 &&
                   csfs[0*norb+2] == 3 && csfs[0*norb+3] == 3);
        checks++;
        printf("  %-55s %s  (csf[0]={%d,%d,%d,%d})\n",
               "csfs[0] == {E,E,F,F} = {0,0,3,3}", ok0 ? "PASS" : "FAIL",
               csfs[0], csfs[1], csfs[2], csfs[3]);
        if (!ok0) failures++;

        /* CSF 1: e,u,d,f */
        int ok1 = (csfs[1*norb+0] == 0 && csfs[1*norb+1] == 1 &&
                   csfs[1*norb+2] == 2 && csfs[1*norb+3] == 3);
        checks++;
        printf("  %-55s %s  (csf[1]={%d,%d,%d,%d})\n",
               "csfs[1] == {E,U,D,F} = {0,1,2,3}", ok1 ? "PASS" : "FAIL",
               csfs[norb], csfs[norb+1], csfs[norb+2], csfs[norb+3]);
        if (!ok1) failures++;
    }

    /* --- hexpr_drt_node_info --- */
    {
        int orb, nve, is;
        hexpr_status_t st = hexpr_drt_node_info(drt, 0, &orb, &nve, &is);
        check_ok("hexpr_drt_node_info(0) returns OK", st);
        check_int("node 0 orb == 0", orb, 0);
        check_int("node 0 nve == 0", nve, 0);
        check_int("node 0 is  == 0", is,  0);

        /* out-of-range: expect HEXPR_ERR_RANGE */
        hexpr_status_t bad = hexpr_drt_node_info(drt, 9999, &orb, &nve, &is);
        checks++;
        int ok = (bad == HEXPR_ERR_RANGE);
        printf("  %-55s %s  (status %d)\n",
               "hexpr_drt_node_info(9999) == HEXPR_ERR_RANGE", ok ? "PASS" : "FAIL", (int)bad);
        if (!ok) failures++;
    }

    /* --- hexpr_expr_build --- */
    hexpr_expr_t *expr = hexpr_expr_build(drt);
    check_notnull("hexpr_expr_build not NULL", expr);
    if (!expr) { hexpr_drt_destroy(drt); return 1; }

    /* --- nterms_full / _one: cross-check against show output --- */
    int64_t nfull = hexpr_expr_nterms_full(expr);
    int64_t none  = hexpr_expr_nterms_one(expr);
    {
        checks++;
        int ok = (nfull > 0);
        printf("  %-55s %s  (nterms_full=%" PRId64 ")\n",
               "nterms_full > 0", ok ? "PASS" : "FAIL", nfull);
        if (!ok) failures++;
    }
    {
        checks++;
        int ok = (none > 0 && none <= nfull);
        printf("  %-55s %s  (nterms_one=%" PRId64 ")\n",
               "0 < nterms_one <= nterms_full", ok ? "PASS" : "FAIL", none);
        if (!ok) failures++;
    }

    /* --- hexpr_expr_read FULL --- */
    int32_t *ij   = malloc((size_t)nfull * sizeof(int32_t));
    int32_t *pqrs = malloc((size_t)nfull * sizeof(int32_t));
    double  *coef = malloc((size_t)nfull * sizeof(double));
    check_ok("hexpr_expr_read(FULL) returns OK",
             hexpr_expr_read(expr, HEXPR_EXPR_FULL, ij, pqrs, coef));

    /* Count 1-electron terms manually and compare to nterms_one */
    int norb = hexpr_drt_norb(drt);
    int n_oel = norb * (norb + 1) / 2;
    int64_t counted_one = 0;
    for (int64_t i = 0; i < nfull; i++)
        if (pqrs[i] < n_oel) counted_one++;
    {
        checks++;
        int ok = (counted_one == none);
        printf("  %-55s %s  (counted=%" PRId64 " nterms_one=%" PRId64 ")\n",
               "1-electron count matches nterms_one", ok ? "PASS" : "FAIL",
               counted_one, none);
        if (!ok) failures++;
    }

    /* --- hexpr_expr_read ONE --- */
    int32_t *ij1   = malloc((size_t)none * sizeof(int32_t));
    int32_t *pqrs1 = malloc((size_t)none * sizeof(int32_t));
    double  *coef1 = malloc((size_t)none * sizeof(double));
    check_ok("hexpr_expr_read(ONE) returns OK",
             hexpr_expr_read(expr, HEXPR_EXPR_ONE, ij1, pqrs1, coef1));
    /* All returned pqrs must be < n_oel */
    {
        int all_ok = 1;
        for (int64_t i = 0; i < none; i++)
            if (pqrs1[i] >= n_oel) { all_ok = 0; break; }
        checks++;
        printf("  %-55s %s\n",
               "hexpr_expr_read(ONE): all pqrs < n_oel", all_ok ? "PASS" : "FAIL");
        if (!all_ok) failures++;
    }

    /* --- Cross-check show_full output against hexpr_expr_read ---
     * Write show_full to tmp file and compare first/last numeric values. */
    system("mkdir -p " TMP_DIR);
    char show_path[512];
    snprintf(show_path, sizeof(show_path), "%s/s_full_4o4e.expr.txt", TMP_DIR);
    FILE *f = fopen(show_path, "w");
    if (f) {
        hexpr_expr_show_full(expr, f);
        fclose(f);
        /* read first term from show output and compare to expr_read */
        /* Reference: first data line: "       1             5      2.000... " -> ij=1, pqrs=5, coef~=2 */
        int mm = compare_files(REF_EXPR, show_path, 1e-10, 1e-10, 3);
        checks++;
        printf("  %-55s %s\n",
               "show_full output matches reference", mm == 0 ? "PASS" : "FAIL");
        if (mm) failures++;
    }

    /* Check first term from expr_read matches reference (ij=0 stored, +1 for display) */
    /* Reference first term: IJ=1 -> ij stored = 0, pqrs=5, coef~=2.0 */
    if (nfull > 0) {
        checks++;
        int ok = (ij[0] == 0 && pqrs[0] == 5 && fabs(coef[0] - 2.0) < 1e-9);
        printf("  %-55s %s  (ij=%d pqrs=%d coef=%.4f)\n",
               "first term: ij=0,pqrs=5,coef~=2.0", ok ? "PASS" : "FAIL",
               ij[0], pqrs[0], coef[0]);
        if (!ok) failures++;
    }

    free(ij); free(pqrs); free(coef);
    free(ij1); free(pqrs1); free(coef1);

    /* --- Cleanup --- */
    hexpr_expr_destroy(expr);
    hexpr_drt_destroy(drt);

    printf("%s (%d/%d checks, %d failure(s))\n",
           failures ? "FAIL" : "PASS", checks - failures, checks, failures);
    return failures ? 1 : 0;
}
