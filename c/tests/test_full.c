/* test_full.c -- end-to-end regression test.
 *
 * Runs run_reference into ./tmp/test_full/ (inside the repo, not /tmp),
 * then compares all 36 generated files against reference/.
 * Cleans up on success; leaves directory on failure for inspection.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "util/compare.h"

#define REF_DIR  "../../reference"
#define TMP_DIR  "tmp/test_full"
#define RUN_REF  "../run_reference"
#define LIB_DIR  ".."
#ifdef __APPLE__
#define HEXPR_LIBENV "DYLD_LIBRARY_PATH="
#else
#define HEXPR_LIBENV "LD_LIBRARY_PATH="
#endif
#define ATOL     1e-10
#define RTOL     1e-10

static const char *CASE_NAMES[] = {
    "s_full_4o4e", "t_full_4o4e",
    "s_full_6o6e", "t_full_6o6e",
    "s_soci_small", "s_foci_small",
    "t_soci_small", "t_foci_small",
    "s_foci_mid",  "s_soci_mid",
    "s_soci_med",  "s_foci_med",
    "t_soci_med",  "t_foci_med",
};
static const int NCASES = 14;

int main(void)
{
    /* Create output directory */
    system("mkdir -p " TMP_DIR);

    /* Run the driver; redirect its stderr progress to /dev/null for clean output */
    char cmd[1024];
    snprintf(cmd, sizeof(cmd),
             HEXPR_LIBENV LIB_DIR " " RUN_REF " -o " TMP_DIR " 2>/dev/null");
    int rc = system(cmd);
    if (rc != 0) {
        fprintf(stderr, "run_reference exited with code %d\n", rc);
        return 1;
    }

    printf("=== Full end-to-end regression test ===\n");

    int failures = 0;

    for (int i = 0; i < NCASES; i++) {
        const char *name = CASE_NAMES[i];
        char out[512], ref[512];
        int mm;

        /* drt.txt */
        snprintf(out, sizeof(out), "%s/%s.drt.txt",      TMP_DIR, name);
        snprintf(ref, sizeof(ref), "%s/%s.drt.txt",      REF_DIR, name);
        if (access(ref, R_OK) != 0) {
            printf("  %-20s  drt      SKIP (reference file not bundled)\n", name);
        } else {
            mm = compare_files(ref, out, ATOL, RTOL, 3);
            printf("  %-20s  drt      %s", name, mm == 0 ? "PASS" : "FAIL");
            if (mm) { printf(" (%d mismatch(es))", mm); failures++; }
            printf("\n");
        }

        /* expr.txt */
        snprintf(out, sizeof(out), "%s/%s.expr.txt",     TMP_DIR, name);
        snprintf(ref, sizeof(ref), "%s/%s.expr.txt",     REF_DIR, name);
        if (access(ref, R_OK) != 0) {
            printf("  %-20s  expr     SKIP (reference file not bundled)\n", name);
        } else {
            mm = compare_files(ref, out, ATOL, RTOL, 3);
            printf("  %-20s  expr     %s", name, mm == 0 ? "PASS" : "FAIL");
            if (mm) { printf(" (%d mismatch(es))", mm); failures++; }
            printf("\n");
        }

        /* expr_one.txt */
        snprintf(out, sizeof(out), "%s/%s.expr_one.txt", TMP_DIR, name);
        snprintf(ref, sizeof(ref), "%s/%s.expr_one.txt", REF_DIR, name);
        if (access(ref, R_OK) != 0) {
            printf("  %-20s  expr_one SKIP (reference file not bundled)\n", name);
        } else {
            mm = compare_files(ref, out, ATOL, RTOL, 3);
            printf("  %-20s  expr_one %s", name, mm == 0 ? "PASS" : "FAIL");
            if (mm) { printf(" (%d mismatch(es))", mm); failures++; }
            printf("\n");
        }
    }

    printf("%s (%d failure(s))\n", failures ? "FAIL" : "PASS", failures);

    if (failures == 0)
        system("rm -rf " TMP_DIR);

    return failures ? 1 : 0;
}
