/* test_expr.c -- compare expr.txt and expr_one.txt against reference */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "../include/hexpr.h"
#include "util/compare.h"

#define REF_DIR  "../../reference"
#define TMP_DIR  "tmp/test_expr"
#define ATOL 1e-10
#define RTOL 1e-10

typedef struct {
    const char *name;
    int  nve, spin, ncore, nval, next;
    hexpr_drt_kind_t kind;
} tcase_t;

static const tcase_t CASES[] = {
    {"s_full_4o4e",  4, 0, 0, 4, 0, HEXPR_DRT_SOCI},
    {"t_full_4o4e",  4, 2, 0, 4, 0, HEXPR_DRT_SOCI},
    {"s_full_6o6e",  6, 0, 0, 6, 0, HEXPR_DRT_SOCI},
    {"t_full_6o6e",  6, 2, 0, 6, 0, HEXPR_DRT_SOCI},
    {"s_soci_small", 6, 0, 1, 4, 1, HEXPR_DRT_SOCI},
    {"s_foci_small", 6, 0, 1, 4, 1, HEXPR_DRT_FOCI},
    {"t_soci_small", 6, 2, 1, 4, 1, HEXPR_DRT_SOCI},
    {"t_foci_small", 6, 2, 1, 4, 1, HEXPR_DRT_FOCI},
    {"s_soci_med",  10, 0, 2, 4, 4, HEXPR_DRT_SOCI},
    {"s_foci_med",  10, 0, 2, 4, 4, HEXPR_DRT_FOCI},
    {"t_soci_med",  10, 2, 2, 4, 4, HEXPR_DRT_SOCI},
    {"t_foci_med",  10, 2, 2, 4, 4, HEXPR_DRT_FOCI},
};
static const int NCASES = (int)(sizeof(CASES)/sizeof(CASES[0]));

static int run_case(const tcase_t *tc)
{
    hexpr_drt_t *drt = hexpr_drt_create(tc->kind, tc->nve, tc->spin,
                                        tc->ncore, tc->nval, tc->next);
    if (!drt) { fprintf(stderr, "[%s] drt_create failed\n", tc->name); return 2; }

    hexpr_expr_t *expr = hexpr_expr_build(drt);
    hexpr_drt_destroy(drt);
    if (!expr) { fprintf(stderr, "[%s] expr_build failed\n", tc->name); return 2; }

    int fail = 0;

    /* expr.txt */
    char out_path[512], ref_path[512];
    snprintf(out_path, sizeof(out_path), "%s/%s.expr.txt", TMP_DIR, tc->name);
    snprintf(ref_path, sizeof(ref_path), "%s/%s.expr.txt", REF_DIR, tc->name);
    FILE *f = fopen(out_path, "w");
    if (!f) { perror(out_path); hexpr_expr_destroy(expr); return 2; }
    hexpr_expr_show_full(expr, f);
    fclose(f);
    if (access(ref_path, R_OK) != 0) {
        printf("  %-20s  expr     SKIP (reference file not bundled)\n", tc->name);
    } else {
        int mm = compare_files(ref_path, out_path, ATOL, RTOL, 3);
        printf("  %-20s  expr     %s", tc->name, mm == 0 ? "PASS" : "FAIL");
        if (mm) printf(" (%d mismatch(es))", mm);
        printf("\n");
        if (mm) fail++;
    }

    /* expr_one.txt */
    snprintf(out_path, sizeof(out_path), "%s/%s.expr_one.txt", TMP_DIR, tc->name);
    snprintf(ref_path, sizeof(ref_path), "%s/%s.expr_one.txt", REF_DIR, tc->name);
    f = fopen(out_path, "w");
    if (!f) { perror(out_path); hexpr_expr_destroy(expr); return 2; }
    hexpr_expr_show_one(expr, f);
    fclose(f);
    if (access(ref_path, R_OK) != 0) {
        printf("  %-20s  expr_one SKIP (reference file not bundled)\n", tc->name);
    } else {
        int mm = compare_files(ref_path, out_path, ATOL, RTOL, 3);
        printf("  %-20s  expr_one %s", tc->name, mm == 0 ? "PASS" : "FAIL");
        if (mm) printf(" (%d mismatch(es))", mm);
        printf("\n");
        if (mm) fail++;
    }

    hexpr_expr_destroy(expr);
    return fail;
}

int main(void)
{
    hexpr_init();
    system("mkdir -p " TMP_DIR);
    printf("=== Expression regression test ===\n");
    int failures = 0;
    for (int i = 0; i < NCASES; i++)
        failures += run_case(&CASES[i]);
    printf("%s (%d failure(s))\n", failures ? "FAIL" : "PASS", failures);
    return failures ? 1 : 0;
}
